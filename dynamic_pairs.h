using namespace psi;   
class sort_pred;
void construct_fock(double *FA,SharedMatrix density,int i,int j,int k,int l,double val,int nao);
void dynamic_pairs(boost::shared_ptr<IntegralFactory> integral,boost::shared_ptr<BasisSet> aoBasis,double *FA,SharedMatrix density,int rank,int numtasks)
{     
        int i,j,k,l,tag=0;
        int nshells = aoBasis->nshell();
        int nbf[] = {aoBasis->nbf()};
        int nao = nbf[0];
        int Nshells = nshells *(nshells + 1)/2 ;
        boost::shared_ptr<TwoBodyAOInt> eri(integral->eri());
        const double *buffer = eri->buffer();

        if (rank==0)
        {
    int finished = 0;
    std::vector<std::pair<int, int> > vec;
    for(i=0; i < nshells; i++)
        for(j=0; j <= i; j++)
        {
            vec.push_back(std::make_pair(i,j));
        }

    std::sort(vec.begin(), vec.end(), sort_pred(aoBasis));

    /*for(int i=0;i<vec.size();i++)
    {
       int a = aoBasis->shell(vec[i].first).nfunction() ;
       int b = aoBasis->shell(vec[i].second).nfunction() ;
    }
    */



    int task_index =0;
    int remaining_requests = Nshells + numtasks -1 ;

        while(remaining_requests)
    {

        int message;
        int pid = (remaining_requests % (numtasks -1)) + 1 ;
        MPI_Recv(&message,1,MPI_INT,pid,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        //printf("\nrank %d received message: %d from rank %d\n",rank,message,pid);

        if (task_index < Nshells)
        {
            int num[3] = {vec[task_index].first,vec[task_index].second,finished};
            MPI_Send(&num[0],3, MPI_INT,pid, tag, MPI_COMM_WORLD);
            task_index +=1 ;
        }

        else
        {
        finished = 1;
        int num[3] = {0,0,finished};
        MPI_Send(&num[0],3, MPI_INT,pid, tag, MPI_COMM_WORLD);
        //printf("\nrank %d sent finished value: %d to rank %d \n",rank,finished,pid);
        }
        remaining_requests = remaining_requests - 1 ;
        }
    }
    
    else
        {
            int i,j;
            int finished =0 ;
            int num[3];
            int message = 1;
            MPI_Send(&message,1, MPI_INT,0, tag, MPI_COMM_WORLD);
            //printf("\nrank %d sent message: %d to rank 0 \n",rank,message);
            MPI_Recv(&num[0],3,MPI_INT,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            i = num[0];
            j = num[1];
            finished = num[2];
            //printf("\nrank %d received finished value: %d from rank 0 \n",rank,finished);
            while(!finished)
            {
                for(k=0; k <= i; k++)
                for(l=0; l <= (i==k ? j : k); l++) {
                eri->compute_shell(i,j,k,l);
                AOIntegralsIterator intIter = integral->integrals_iterator(i,j,k,l);
                for(intIter.first(); intIter.is_done() == false; intIter.next())
                {
                    int i = intIter.i();
                    int j = intIter.j();
                    int k = intIter.k();
                    int l = intIter.l();
                    double val  = buffer[intIter.index()];
                    construct_fock(FA,density,i,j,k,l,val,nao);

                }
            }

            MPI_Send(&message,1, MPI_INT,0, tag, MPI_COMM_WORLD);
            //printf("\nrank %d sent message: %d to rank 0 \n",rank,message);
            MPI_Recv(&num[0],3,MPI_INT,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            i=num[0];
            j=num[1];
            finished = num[2];
            //printf("\nrank %d received finished value: %d from rank 0 \n",rank,finished);
        }
    }
}
