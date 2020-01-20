// ***************************************************************************************************************
//
//          Mini-Aevol is a reduced version of Aevol -- An in silico experimental evolution platform
//
// ***************************************************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <https://gitlab.inria.fr/rouzaudc/mini-aevol>
// Web: https://gitlab.inria.fr/rouzaudc/mini-aevol
// E-mail: See <jonathan.rouzaud-cornabas@inria.fr>
// Original Authors : Jonathan Rouzaud-Cornabas
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ***************************************************************************************************************



#include <iostream>
#include <getopt.h>
#include <cstring>

#include "ExpManager.h"
#include <mpi.h>

void print_help(char* prog_path) {
    // Get the program file-name in prog_name (strip prog_path of the path)
    char* prog_name; // No new, it will point to somewhere inside prog_path
    if ((prog_name = strrchr(prog_path, '/'))) prog_name++;
    else prog_name = prog_path;

    printf("******************************************************************************\n");
    printf("*                                                                            *\n");
    printf("*                   mini-aevol - Artificial Evolution                        *\n");
    printf("*                                                                            *\n");
    printf("* Aevol is a simulation platform that allows one to let populations of       *\n");
    printf("* digital organisms evolve in different conditions and study experimentally  *\n");
    printf("* the mechanisms responsible for the structuration of the genome and the     *\n");
    printf("* transcriptome.                                                             *\n");
    printf("*                                                                            *\n");
    printf("* This is a mini-version with a cutdown biological model and a lot of        *\n");
    printf("* features missing.                                                          *\n");

    printf("*        IF YOU WANT TO USE IT TO STUDY EVOLUTION, DO NOT DO IT !            *\n");
    printf("*                                                                            *\n");
    printf("******************************************************************************\n");
    printf("\n");
    printf("Usage : %s -e or --help\n", prog_name);
    printf("   or : %s \n", prog_name);
    printf("   or : %s -r GENERATION_STEP\n", prog_name);
    printf("\nOptions\n");
    printf("  -e, --help\tprint this help, then exit\n");
    printf("  -n, --nbsteps\tNumber of generation to run\n");
    printf("  -w, --width WIDTH_SIZE\tWidth of the population grid is WIDTH_SIZE\n");
    printf("  -h, --height HEIGHT_SIZE\tHeight of the population grid is HEIGHT_SIZE\n");
    printf("  -m, --mutation_rate MUTATION_RATE\tMutation rate is set to MUTATION_RATE\n");
    printf("  -g, --genome_size GENOME_SIZE\tGenome at the initial genome is GENOME_SIZE bps\n");
    printf("  -b, --backup_step BACKUP_STEP\tDo a simulation backup/checkpoint every BACKUP_STEP\n");
    printf("  -r, --resume RESUME_STEP\tResume the simulation from the RESUME_STEP generations\n");
}

int main(int argc, char* argv[]) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Find out rank, size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // We are assuming at least 2 processes for this task
    if (world_size < 2) {
        fprintf(stderr, "World size must be greater than 1 for %s\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    printf("Processor %d:%s has started\n",world_rank, processor_name);
/*
    int nbstep = -1;
    int width = -1;
    int height = -1;
    double mutation_rate = -1;
    int genome_size = -1;
    int resume = -1;
    int backup_step = -1;
    int seed = -1;
    int nb_threads = -1;

    const char * options_list = "e:::n:w:h:m:g:b:r:s:t:";
    static struct option long_options_list[] = {
            // Print help
            { "help",     no_argument,        NULL, 'e' },
            // Number of generations to be run
            { "nsteps",  required_argument,  NULL, 'n' },
            // Width size of the grid
            { "width", required_argument,  NULL, 'w' },
            // Height size of the grid
            { "height", required_argument,  NULL, 'h' },
            // Mutation rate
            { "mutation_rate", required_argument,  NULL, 'm' },
            // Size of the initial genome
            { "genome_size", required_argument,  NULL, 'g' },
            // Resuming from generation X
            { "resume", required_argument,  NULL, 'r' },
            // Backup step
            { "backup_step", required_argument,  NULL, 'b' },
            // Seed
            { "seed", required_argument,  NULL, 's' },
            { "nb_threads", required_argument,  NULL, 't' },
            { 0, 0, 0, 0 }
    };


    // -------------------------------------------------------------------------
    // 3) Get actual values of the command-line options
    // -------------------------------------------------------------------------
    int option;
    while ((option =
                    getopt_long(argc, argv, options_list, long_options_list, NULL))
           != -1) {
        switch (option) {
            case 'e' : {
                print_help(argv[0]);
                exit(EXIT_SUCCESS);
            }
            case 'w' : {
                width = atoi(optarg);
                break;
            }
            case 'h' : {
                height = atoi(optarg);
                break;
            }
            case 'm' : {
                mutation_rate = atof(optarg);
                break;
            }
            case 'g' : {
                genome_size = atoi(optarg);
                break;
            }
            case 'r' : {
                resume = atoi(optarg);
                break;
            }
            case 'b' : {
                backup_step = atoi(optarg);
                break;
            }
            case 's' : {
                seed = atoi(optarg);
                break;
            }
            case 'n' : {
                nbstep = atoi(optarg);
                break;
            }
            case 't' : {
                nb_threads = atoi(optarg);
                break;
            }
            default : {
                // An error message is printed in getopt_long, we just need to exit
                printf("Error unknown parameter\n");
                exit(EXIT_FAILURE);
            }
        }
    }


    if (resume >= 0) {
        if ((width != -1) || (height != -1)|| (mutation_rate != -1.0) || (genome_size != -1) ||
            (backup_step != -1) || (seed != -1)) {
            printf("Parameter(s) can not change during the simulation (i.e. when resuming a simulation, parameter(s) can not change)\n");
            exit(EXIT_FAILURE);
        }
        if (nbstep == -1) nbstep = 1000;
    } else {
        if (nbstep == -1) nbstep = 1000;
        if (width == -1) width = 32;
        if (height == -1) height = 32;
        if (mutation_rate == -1) mutation_rate = 0.00001;
        if (genome_size == -1) genome_size = 5000;
        if (backup_step == -1) backup_step = 1000;
        if (seed == -1) seed = 566545665;
        if (nb_threads == -1) nb_threads = 1;
    }

*/

std::ofstream res_file;
if (world_rank == 0){
    res_file.open("./simulation_example_/results_nb_gens.csv", std::ios_base::app);
}

for (int nb_threads = 1; nb_threads <= 4; nb_threads++){
    for (int nb_gens = 1000; nb_gens <= 20000; nb_gens += 1000){
        ExpManager *exp_manager;
        exp_manager = new ExpManager(32,32, 1337, 0.00001, 4096, 0.03, 1000, 1000, nb_threads);
        auto t1 = std::chrono::high_resolution_clock::now();
        exp_manager->run_evolution(nb_gens);
        auto t2 = std::chrono::high_resolution_clock::now();
        if(world_rank == 0 && nb_gens % 1000 == 0){
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
            res_file << 'mpi,' << 2 * nb_threads << ',' << nb_gens << ',' << duration << '\n';
        }
        delete exp_manager;
    }

    /*
    for (int genom_length = 1024; genom_length <= 8192; genom_length = genom_length * 2){
        ExpManager *exp_manager;
        exp_manager = new ExpManager(32,32, 1337, 0.00001, genom_length, 0.03, 1000, 1000, nb_threads);
        auto t1 = std::chrono::high_resolution_clock::now();
        exp_manager->run_evolution(10000);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        res_file << nb_threads << ',' << genom_length << ',' << duration << '\n';
        delete exp_manager;
    }
    */

    /* 
    for (double muration_rate = 0.05; muration_rate <= 0.5; muration_rate += 0.05){
        ExpManager *exp_manager;
        exp_manager = new ExpManager(32,32, 1337, muration_rate, 2048, 0.03, 1000, 1000, nb_threads);
        auto t1 = std::chrono::high_resolution_clock::now();
        exp_manager->run_evolution(20);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        res_file << nb_threads << ',' << muration_rate << ',' << duration << '\n';
        delete exp_manager;
    }
    */
    
}

    MPI_Finalize();
    return 0;
}
