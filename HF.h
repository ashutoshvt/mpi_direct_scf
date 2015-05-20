using namespace psi;
void static_pairs(boost::shared_ptr<IntegralFactory>integral,boost::shared_ptr<BasisSet> aoBasis, double *FA,SharedMatrix density,int rank,int numtasks);
void dynamic_pairs(boost::shared_ptr<IntegralFactory>integral,boost::shared_ptr<BasisSet> aoBasis, double *FA,SharedMatrix density,int rank,int numtasks);
void static_quartets(boost::shared_ptr<IntegralFactory> integral,boost::shared_ptr<BasisSet> aoBasis,double *FA, SharedMatrix density,int rank,int numtasks);
class HF
{
 public:

 int nao;
 int rank;
 int numtasks; 
 int natom;
 int occ_double;

 double elec_energy;
 double rmsd_hf;
 Options options;
 std::string guess;
 std::string structure;
 std::string distribution;
 boost::shared_ptr<IntegralFactory> integral;
 boost::shared_ptr<BasisSet> aoBasis;
 boost::shared_ptr<MatrixFactory> factory;
 boost::shared_ptr<Molecule> molecule;

 SharedMatrix sMat;
 SharedMatrix basisno;
 SharedMatrix tMat;
 SharedMatrix vMat;
 SharedMatrix hMat;   
 SharedMatrix fock;   
 SharedMatrix fock_on;   
 SharedMatrix density;   
 SharedMatrix density_old;   
 SharedMatrix coeff;   
 SharedMatrix evecs;   
 SharedVector mo_energy;
 SharedMatrix symm;
 SharedVector evals;

HF( boost::shared_ptr<IntegralFactory> integral,boost::shared_ptr<BasisSet> aoBasis,boost::shared_ptr<MatrixFactory> factory,boost::shared_ptr<Molecule> molecule,Options &options,int rank,int numtasks)  
{

this->aoBasis = aoBasis;
this->integral = integral;
this->factory = factory;
this->molecule = molecule;
this->rank = rank;
this->numtasks = numtasks;
this->options = options;
int nbf[] = {aoBasis->nbf()};
nao = nbf[0];
natom = molecule->natom();
for(int i=0;i< natom;i++)
occ_double += molecule->Z(i);
occ_double /= 2;

guess = options.get_str("GUESS");
distribution = options.get_str("DISTRIBUTION");
structure = options.get_str("STRUCTURE");



}

void allocate_memory()
{
      sMat = SharedMatrix(factory->create_matrix("Overlap"));
      tMat= SharedMatrix(factory->create_matrix("Kinetic"));
      vMat= SharedMatrix(factory->create_matrix("Potential"));
      hMat= SharedMatrix(factory->create_matrix("One Electron Ints"));
      evecs= SharedMatrix(factory->create_matrix(" eigen vectors "));
      symm= SharedMatrix(factory->create_matrix(" symmetric matrix "));
      coeff= SharedMatrix(factory->create_matrix(" mo coefficients matrix "));
      fock= SharedMatrix(factory->create_matrix(" fock matrix "));
      fock_on= SharedMatrix(factory->create_matrix(" orthogonalized fock matrix "));
      basisno= SharedMatrix(factory->create_matrix(" non-orthogonal original basis "));
      density= SharedMatrix(factory->create_matrix("density matrix "));
      density_old= SharedMatrix(factory->create_matrix(" old density matrix "));
      evals = SharedVector(new Vector(" eigen values ",nao));
      mo_energy = SharedVector(new Vector("mo energy",nao));
}


void compute_overlap()
{
   boost::shared_ptr<OneBodyAOInt> sOBI(integral->ao_overlap());
   sOBI->compute(sMat);  
}

void compute_hcore()
{
   boost::shared_ptr<OneBodyAOInt> tOBI(integral->ao_kinetic());
   boost::shared_ptr<OneBodyAOInt> vOBI(integral->ao_potential());
   tOBI->compute(tMat);
   vOBI->compute(vMat);
   hMat->copy(tMat); 
   hMat->add(vMat); 
}


void compute_density()
{
        coeff->zero();
        mo_energy->zero();
        fock_on->zero();
        basisno->zero();
        density_old->zero();
        symm->zero();

        sMat->diagonalize(evecs,evals);

        for( int i=0;i<nao;i++)
        for( int j=0;j<nao;j++)
        for(int k=0;k<nao;k++)
        symm->add(i,j, evecs->get(i,k)*(1.0/sqrt(evals->get(k))) * evecs->get(j,k))  ;

        fock_on->transform(fock,symm);
        fock_on->diagonalize(coeff,mo_energy);
        basisno->gemm(0,0,1,symm,coeff,0);
        density_old->copy(density);
        density->zero();

        for( int i=0;i<nao;i++)
        for(int j=0;j<nao;j++)
        for(int m=0;m< occ_double;m++)
        density->add(i,j,basisno->get(i,m) * basisno->get(j,m));
}



void guess_density()
{
   if (guess == "CORE")
        {
         fock->copy(hMat);
         compute_density();
        } 
    else if (guess == "SAD")
        {
          boost::shared_ptr<scf::SADGuess> guess(new scf::SADGuess(aoBasis,occ_double,occ_double,options));
          guess->compute_guess();
          density->copy(guess->Da());
          fock->copy(hMat);
          compute_density();  
        }
}

double energy()
{
     elec_energy = 0;
    for( int i=0;i<nao;i++)
    for(int j=0;j<nao;j++)
    elec_energy += density->get(i,j) *( hMat->get(i,j) + fock->get(i,j)); 
    return elec_energy;
}

double rmsd()
{
    rmsd_hf = 0;
    for( int i=0;i<nao;i++)
    for(int j=0;j<nao;j++)
    rmsd_hf += (density->get(i,j)- density_old->get(i,j))* (density->get(i,j)- density_old->get(i,j));
    rmsd_hf = sqrt(rmsd_hf)/nao ;
    return rmsd_hf ;
}


void compute_fock()
{

 #define FA(i,j) FA[i*nao+j]

double *FA = new double[nao*nao];

    for (int i=0;i<nao;i++)
    for (int j=0;j<nao;j++)
    FA(i,j) = 0;


    if (structure == "PAIRS" && distribution == "STATIC")
    { static_pairs(integral,aoBasis,FA,density,rank,numtasks);} // pairs, static if ends

    else if (structure == "PAIRS" && distribution == "DYNAMIC")
       { dynamic_pairs(integral,aoBasis,FA,density,rank,numtasks) ;}

    else if (structure == "QUARTET")
    { static_quartets(integral,aoBasis,FA,density,rank,numtasks);}

    MPI_Allreduce(MPI_IN_PLACE, FA, nao*nao, MPI::DOUBLE, MPI::SUM,MPI_COMM_WORLD);


    for (int i=0;i<nao;i++)
    for (int j=0;j<nao;j++)
    fock->set(i,j,FA(i,j));

    fock->add(hMat);



}

};
