using namespace psi;
void static_quartets(boost::shared_ptr<IntegralFactory> integral,boost::shared_ptr<BasisSet> aoBasis,double *FA, SharedMatrix density,int rank,int numtasks)
{
    int nbf[] = { aoBasis->nbf() };
    int nao = nbf[0];
    int count=0;
    boost::shared_ptr<TwoBodyAOInt> eri(integral->eri());
    const double *buffer = eri->buffer();
  AOShellCombinationsIterator shellIter = integral->shells_iterator();
        for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
               if (count % numtasks == rank)
            {
        eri->compute_shell(shellIter);
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next())
         {
                int i = intIter.i();
                int j = intIter.j();
                int k = intIter.k();
                int l = intIter.l();
                double val = buffer[intIter.index()] ;
                construct_fock(FA,density,i,j,k,l,val,nao);
            }
        }
        count ++ ;
    }
}

