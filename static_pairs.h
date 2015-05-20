using namespace psi;
void static_pairs(boost::shared_ptr<IntegralFactory> integral,boost::shared_ptr<BasisSet> aoBasis,double *FA, SharedMatrix density,int rank,int numtasks)
{
        int nbf[] = { aoBasis->nbf() };
        int nao = nbf[0];
        int nshells = aoBasis->nshell();
        boost::shared_ptr<TwoBodyAOInt> eri(integral->eri());
        const double *buffer = eri->buffer();
        int i,j,k,l;
        double val;
        int count =0;
    
    for(i=0; i < nshells; i++)
        for(j=0; j <= i; j++)
        {
          if (count % numtasks == rank)
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
            }
        count ++ ;
    }
}

