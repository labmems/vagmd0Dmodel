/*
Computational code for a simple, reduced model for the analysis of vacuum-enhanced air gap membrane distillation (V-AGMD) desalination modules

Contributors:
Kleber Marques Lisboa, Laboratory of Thermal Sciences (LATERMO), Universidade Federal Fluminense

Acknowledgements:
This work was funded by Petrogal do Brasil S.A. and FAPERJ
*/

static char help[] = "Reduced model for vacuum-enhanced air gap membrane distillation (V-AGMD) modules\n\n"
"Developer: Prof. Kleber Marques Lisboa, UFF\n\n"
"Funding: Petrogal do Brasil S.A. and FAPERJ\n\n"
"Usage:\n\n"
"Let $BINFOLDER be the folder that contains the compiled binary; if running from the project root, $BINFOLDER=./bin/\n"
"In the command line, type $ $BINFOLDER/vagmd0Dmodel -variable1 value1 -variable2 value2\n\n"
"For instance, $ $BINFOLDER/vagmd0Dmodel -feed_mass_flow_rate 0.08 -vacuum_pressure -50000.0 ...\n\n"
"Command-line arguments:\n\n"
"-feed_mass_flow_rate: type double, unit kg/s\n"
"Description - Mass flow rate through the feed channel.\n\n"
"-cool_mass_flow_rate: type double, unit kg/s\n"
"Description - Mass flow rate through the cooling channel\n\n"
"-entry_temperature_feed: type double, unit °C\n"
"Description - Feed temperature at the inlet of the module\n\n"
"-entry_temperature_cool: type double, unit °C\n"
"Description - Coolant temperature at the inlet of the module\n"
"-entry_salinity_feed: type double, unit wt%%\n"
"Description - Feed salinity at the inlet of the module\n\n"
"-entry_salinity_cool: type double, unit wt%%\n"
"Description - Coolant salinity at the inlet of the module\n\n"
"-vacuum_pressure: type double, unit Pa (must be negative!)\n"
"Description - Pressure deduced from the atmospheric within the air gap.\n\n"
"-membrane_area: type double, unit m²\n"
"Description - Area of the membrane.\n\n"
"-membrane_thickness: type double, unit m\n"
"Description - Thickness of the membrane.\n\n"
"-membrane_porosity: type double, unit none\n"
"Description - Porosity of the membrane.\n\n"
"-pore_diameter: type double, unit m\n"
"Description - Average pore diameter of the membrane.\n\n"
"-feed_channel_height: type double, unit m\n"
"Description - Height of the feed channel.\n\n"
"-cold_channel_height: type double, unit m\n"
"Description - Height of the cooling channel.\n\n"
"-channel_width: type double, unit m\n"
"Description - Width of both the feed and cooling channels.\n\n"
"-number_channels: type integer, unit m\n"
"Description - Number of feed (or cooling) channels within the evelope.\n\n"
"-spacer_porosity: type double, unit none\n"
"Description - Porosity the channel's spacer.\n\n"
"-gap_spacer_porosity: type double, unit none\n"
"Description - Porosity of the gap's spacer.\n\n"
"-air_gap_thickness: type double, unit m\n"
"Description - Thickness of the air gap.\n\n"
"-wall_thickness: type double, unit m\n"
"Description - Thickness of the condensing wall.\n\n"
"-polymer_conductivity: type double, unit W/mK\n"
"Description - Thermal conductivity of the polymer from which the membrane is made of.\n\n"
"-spacer_conductivity: type double, unit W/mK\n"
"Description - Thermal conductivity of the material from which the spacer is made of.\n\n"
"-wall_conductivity: type double, unit W/mK\n"
"Description - Thermal conductivity of the condensing wall.\n\n";

#include "lib.h"

int main(int argc, char **argv)
{
    PetscMPIInt size;
    EntryData entry_data;

    //-----------------------------------------------------------------------------------------------------------------------------------------------//
    // Initializing PETSc                                                                                                                            //
    //-----------------------------------------------------------------------------------------------------------------------------------------------//

    PetscFunctionBeginUser;
    PetscInitialize(&argc, &argv, (char *)0, help);
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCheck(size == 1, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE, "This program is intended for serial mode only!\n");

    //-----------------------------------------------------------------------------------------------------------------------------------------------//
    // Fetching the entry data                                                                                                                       //
    //-----------------------------------------------------------------------------------------------------------------------------------------------//

    PetscCall(EntryDataBuild(&entry_data));

    //-----------------------------------------------------------------------------------------------------------------------------------------------//
    // Running the model of the plant                                                                                                                //
    //-----------------------------------------------------------------------------------------------------------------------------------------------//

    PetscCall(RunPlant(&entry_data));

    //-----------------------------------------------------------------------------------------------------------------------------------------------//
    // Finalizing PETSc and the program                                                                                                              //
    //-----------------------------------------------------------------------------------------------------------------------------------------------//

    PetscFinalize();

    return 0;
}