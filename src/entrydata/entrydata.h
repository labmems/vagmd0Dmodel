#ifndef ENTRYDATA

#define ENTRYDATA

#include <petscsnes.h>
#include <petscdm.h>
#include <petscdmda.h>

/*
Defining necessary constants
*/

static const PetscReal water_molar_mass = 18.015e-3;
static const PetscReal salt_molar_mass = 58.443e-3;
static const PetscReal gas_constant = 8.3144698;
static const PetscReal atm_pressure = 101.325e3;
static const PetscReal membrane_tortuosity = 2.27; // https://doi.org/10.1016/j.memsci.2017.04.002

// Data structure containing the data involved in the model for the desalination module
typedef struct
{
    // Operational conditions
    PetscReal feed_mass_flow_rate, cool_mass_flow_rate, entry_temperature_feed, entry_temperature_cool,
              entry_salinity_feed, entry_salinity_cool, vacuum_pressure;

    // Geometrical dimensions and fixed properties
    PetscReal membrane_area, membrane_thickness, membrane_porosity, pore_diameter, feed_channel_height,
              cool_channel_height, channel_width, spacer_porosity, gap_spacer_porosity, air_gap_thickness,
              wall_thickness, polymer_conductivity, spacer_conductivity, wall_conductivity;
    PetscInt number_channels;

    // Iterative data
    PetscReal out_temperature_feed, out_temperature_cool, feed_membrane_temperature, gap_membrane_temperature,
              film_boundary_temperature, film_wall_temperature, cool_wall_temperature, out_salinity_feed,
              mass_flux, heat_flux, vapor_heat_flux, feed_outflow_rate;
} DessalData;


// Aggregate of entry data for all components of the system
typedef struct
{
    DessalData dessal_data;
} EntryData;

// Entry data constructor
PetscErrorCode EntryDataBuild(EntryData *entry_data);

#endif