#ifndef PHYSICS

#define PHYSICS

#include "../properties/properties.h"

// Function to calculate the heat transfer coefficients in the water channels
PetscReal ChannelHeatTransfCoef(SaltWaterProperties *bulk_water_prop,
                                SaltWaterProperties *wall_water_prop,
                                PetscReal mass_flow_rate,
                                PetscReal channel_height,
                                PetscReal channel_width,
                                PetscInt number_channels,
                                PetscReal spacer_porosity);

// Function to calculate the effective thermal conductivity of the membrane
PetscReal MembraneConductivity(MoistAirProperties *pore_air_prop,
                               PetscReal polymer_conductivity,
                               PetscReal membrane_porosity);

// Function to calculate the distillate mass flux across the membrane
PetscReal MassFlux(PetscReal membrane_porosity,
                   PetscReal membrane_tortuosity,
                   PetscReal membrane_thickness,
                   PetscReal pore_diameter,
                   PetscReal gap_spacer_porosity,
                   PetscReal air_gap_thickness,
                   PetscReal temperature_membrane,
                   PetscReal temperature_gap,
                   PetscReal feed_membrane_pressure,
                   PetscReal film_boundary_pressure,
                   PetscReal vacuum_pressure);

#endif