/*
 *@BEGIN LICENSE
 *
 * mpi_direct_scf by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <mpi.h>
#include <libqt/qt.h>
#include <libscf_solver/sad.h>
#include <libscf_solver/hf.h>
#include <libscf_solver/rhf.h>
#include "walltime.h"
#include "sort.h"
#include "fock.h"
#include "HF.h"
#include "static_pairs.h"
#include "static_quartets.h"
#include "dynamic_pairs.h"

INIT_PLUGIN

double wall_time();

namespace psi{ namespace mpi_direct_scf {

extern "C"
int read_options(std::string name, Options &options)
{
    if (name == "MPI_DIRECT_SCF"|| options.read_globals()) {
        options.add_int("PRINT", 1);
        options.add_bool("DO_TEI", true);
        options.add_str("STRUCTURE", "quartet");
        options.add_str("GUESS", "CORE");
        options.add_str("GUESS", "SAD");
        options.add_str("STRUCTURE", "pairs");
        options.add_str("DISTRIBUTION", "static");
        options.add_str("DISTRIBUTION", "dynamic");
        options.add_int("MAXITER",50);
        options.add_double("CONV",1e-13);
    }

    return true;
}

extern "C"
PsiReturnType mpi_direct_scf(Options &options)
{

    int numtasks, len, rank, dest,source, rc ;
    char inmsg, outmsg='x';
    MPI_Status Stat;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(hostname, &len);
    printf ("Number of tasks= %d My rank= %d Running on %s\n", numtasks,rank,hostname);
    printf("Task %d starting...\n",rank);

    int print = options.get_int("PRINT");
    int doTei = options.get_bool("DO_TEI");
    int maxiter = options.get_int("MAXITER");
    double conv = options.get_double("CONV");
    std::string structure = options.get_str("STRUCTURE");
    std::string guess = options.get_str("GUESS");
    std::string distribution = options.get_str("DISTRIBUTION");

    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    boost::shared_ptr<BasisSet> basisset = BasisSet::pyconstruct_orbital(Process::environment.molecule(),
            "BASIS", options.get_str("BASIS"));
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory
            (aoBasis, aoBasis, aoBasis, aoBasis));
    // N.B. This should be called after the basis has been built, because the geometry has not been
    // fully initialized until this time.
    molecule->print();
    int nbf[] = {aoBasis->nbf()} ;
    int nao = nbf[0];
    int iter = 0;
    double nucrep = molecule->nuclear_repulsion_energy();
    psi::outfile->Printf("\n    Nuclear repulsion energy: %16.8f\n\n", nucrep);
    boost::shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(1, nbf, nbf);

    SharedVector energy_elec = SharedVector( new Vector(" Electron energy ",maxiter));
    SharedVector delta_e = SharedVector( new Vector(" change in energy ", maxiter));
    SharedVector rmsd = SharedVector( new Vector(" change in root mean square density",maxiter));

    HF hf(integral,aoBasis,factory,molecule,options,rank,numtasks);
    hf.allocate_memory();
    hf.compute_overlap();
    hf.compute_hcore();
    hf.guess_density();
    energy_elec->set(iter,hf.energy());
        

    if (rank ==0) {printf("electronic energy : %lf\n\n",energy_elec->get(iter) + nucrep) ;}

    for(iter=1;iter<maxiter;iter++)
         {

        hf.compute_fock();
        hf.compute_density();
        energy_elec->set(iter,hf.energy());
        rmsd->set(iter,hf.rmsd()) ;
        delta_e->set(iter,energy_elec->get(iter)- energy_elec->get(iter -1 ));
   
        if (rank==0)
        {psi::outfile->Printf("\n %d %20.12lf %20.12lf %20.12lf ", iter,energy_elec->get(iter)+nucrep,delta_e->get(iter),rmsd->get(iter));}
        bool converged = fabs(delta_e->get(iter)) < conv && fabs(rmsd->get(iter)) < conv ;   
         MPI_Bcast(&converged, sizeof(converged), MPI::BYTE, 0,MPI_COMM_WORLD);
            if (converged) {if (rank == 0) printf("\n Converged!\n");
                            break;}
    }
    return Success;
}

}} // End Namespaces
