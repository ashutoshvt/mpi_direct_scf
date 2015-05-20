void construct_fock(double *FA,psi::SharedMatrix density,int i,int j,int k,int l,double val,int nao)
    {
        
        #define FA(i,j) FA[i*nao+j]
    

        int idl1,idl2,unq1,unq2;

        if((i==j)&&(j==k)&&(k==l))
        {
        FA(i,i)+= density->get(i,i) * val ;
        }
         else if ( ((i==j) &&(k!=l) &&(k!=i)) || ((k==l)&&(i!=j)&&(j!=k))  || ((i==k)&&(j==l)&&(j!=k)) || ( (i==j)&&(j==k)&&(i!=l)) || ((j==k)&&(k==l)&&(j!=i)))
        {
        if (i==j)
        {  idl1 = i;
           idl2 = i;
           unq1 = k;
           unq2 = l;
        }
         else if( k==l)
        {
          idl1 = k;
          idl2 = k;
          unq1 = i;
          unq2 = j;
        }
         else if((i==k) &&(j==l))
        {
        idl1 = i;
        idl2 = j;
        unq1 = i;
        unq2 = j;
        }
              FA(idl1,idl2)+= 2* density->get(unq1,unq2) * val;
              FA(idl2,idl1)+= 2* density->get(unq2,unq1) * val;
              FA(idl1,unq1)+=-(density->get(idl2,unq2)*val);
              FA(idl2,unq2)+= -(density->get(idl1,unq1) * val);
              FA(unq1,unq2)+= 2* density->get(idl2,idl1) * val;
              FA(unq2,unq1)+= 2* density->get(idl1,idl2) * val;
              FA(unq1,idl2)+=-(density->get(unq2,idl1)* val);
              FA(unq2,idl1)+= -(density->get(unq1,idl2) * val);
        }

                    // i==j && k==l

        else if((i==j)&&(k==l))
        {
        FA(i,i)+=2*density->get(k,k)*val;
        FA(k,k)+=2*density->get(i,i)*val;
        FA(i,k)+=-(density->get(i,k)*val);
        FA(k,i)+=-(density->get(k,i)*val);
        }
        
        else
        {
              FA(i,j)+= 2* density->get(k,l) * val;
              FA(i,k)+=-(density->get(j,l) * val);
              // now the different permutationally similar integrals //
              // (ij||lk)
              FA(i,j)+= 2* density->get(l,k) * val;
              FA(i,l)+=-(density->get(j,k)*val);
              // (ji||kl)
              FA(j,i)+= 2* density->get(k,l) * val;
              FA(j,k)+=-(density->get(i,l) * val);
              // (ji||lk)
               FA(j,i)+= 2* density->get(l,k) * val;
              FA(j,l)+=-(density->get(i,k)* val);
              // (kl||ij)
              FA(k,l)+= 2* density->get(i,j) * val;
              FA(k,i)+=-(density->get(l,j)* val);
              // (lk||ij)
              FA(l,k)+= 2* density->get(i,j) * val;
              FA(l,i)+=-(density->get(k,j)*val);
              // (kl||ji)
              FA(k,l)+= 2* density->get(j,i) * val;
              FA(k,j)+=-(density->get(l,i) *val);
              // (lk||ji)
              FA(l,k)+= 2* density->get(j,i) * val;
              FA(l,j)+=-(density->get(k,i)*val);
        }   
    }
