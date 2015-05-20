    
using namespace psi;

void compute_hf_density(SharedMatrix sMat, SharedMatrix fock, SharedMatrix density, boost::shared_ptr<MatrixFactory> factory)
{ 
       
    SharedMatrix evecs(factory->create_matrix(" eigen vectors "));
    SharedMatrix symm(factory->create_matrix(" symmetric matrix "));
    SharedMatrix coeff(factory->create_matrix(" mo coefficients matrix "));
    SharedMatrix fock(factory->create_matrix(" fock matrix "));
    SharedMatrix fock_on(factory->create_matrix(" orthogonalized fock matrix "));
    SharedMatrix basisno(factory->create_matrix(" non-orthogonal original basis "));
    SharedMatrix density(factory->create_matrix("density matrix "));
    SharedMatrix density_old(factory->create_matrix(" old density matrix "));
    SharedMatrix temp(factory->create_matrix(" temporary matrix "));
    SharedVector evals = SharedVector( new Vector(" eigen values ",nao));
    SharedVector mo_energy = SharedVector( new Vector("mo energy",nao));


        temp->zero();
        coeff->zero();
        mo_energy->zero();
        fock_on->zero();
        basisno->zero();
        density_old->zero();

        sMat-> diagonalize( evecs,evals);
        symm->zero();

        for( int i=0;i<nao;i++)
        for( int j=0;j<nao;j++)
        for(int k=0;k<nao;k++)
        symm-> add(i,j, evecs->get(i,k)*(1.0/sqrt(evals->get(k))) * evecs->get(j,k))  ;

        temp->gemm(1,0,1,symm,fock,0);
        fock_on->gemm(0,0,1,temp,symm,0);
        fock_on->diagonalize(coeff,mo_energy);
        basisno->gemm(0,0,1,symm,coeff,0);
        density_old->copy(density);
        density->zero();

        for( int i=0;i<nao;i++)
        for(int j=0;j<nao;j++)
        for(int m=0;m< occ_double;m++)
        density->add(i,j,basisno->get(i,m) * basisno->get(j,m));
}
