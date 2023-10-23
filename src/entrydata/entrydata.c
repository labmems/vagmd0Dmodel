#include "entrydata.h"

PetscErrorCode EntryDataBuild(EntryData *entry_data)
{
    PetscFunctionBeginUser;

    DessalData dessal_data;
    // Operational conditions
    PetscReal feed_mass_flow_rate = 400.0 / 3600.0, // Default: 400 kg/h
              cool_mass_flow_rate = 400.0 / 3600.0, // Default: 400 kg/h
              entry_temperature_feed = 60.0, // Default: 60 degC
              entry_temperature_cool = 25.0, // Default: 25.0 degC
              entry_salinity_feed = 3.5e-2, // Default: 3.5 wt%
              entry_salinity_cool = 3.5e-2, // Default: 3.5 wt%
              vacuum_pressure = -50000.0; // Default: -50000.0 Pa

    PetscOptionsGetReal(NULL, NULL, "-feed_mass_flow_rate", &feed_mass_flow_rate, NULL);
    PetscOptionsGetReal(NULL, NULL, "-cool_mass_flow_rate", &cool_mass_flow_rate, NULL);
    PetscOptionsGetReal(NULL, NULL, "-entry_temperature_feed", &entry_temperature_feed, NULL);
    PetscOptionsGetReal(NULL, NULL, "-entry_temperature_cool", &entry_temperature_cool, NULL);
    PetscOptionsGetReal(NULL, NULL, "-entry_salinity_feed", &entry_salinity_feed, NULL);
    PetscOptionsGetReal(NULL, NULL, "-entry_salinity_cool", &entry_salinity_cool, NULL);
    PetscOptionsGetReal(NULL, NULL, "-vacuum_pressure", &vacuum_pressure, NULL);

    dessal_data.feed_mass_flow_rate = feed_mass_flow_rate;
    dessal_data.cool_mass_flow_rate = cool_mass_flow_rate;
    dessal_data.entry_temperature_feed = entry_temperature_feed;
    dessal_data.entry_temperature_cool = entry_temperature_cool;
    dessal_data.entry_salinity_feed = entry_salinity_feed;
    dessal_data.entry_salinity_cool = entry_salinity_cool;
    dessal_data.vacuum_pressure = vacuum_pressure;

    // Geometrical dimensions and fixed properties
    PetscReal membrane_area = 12.96, // Default: 12.96 m^2
              membrane_thickness = 100.0e-6, // Default: 100 microns
              membrane_porosity = 0.85, // Default: 85%
              pore_diameter = 0.32e-6, // Default: 0.32 microns
              feed_channel_height = 2.0e-3, // Default: 2 mm
              cool_channel_height = 2.0e-3, // Default: 2 mm
              channel_width = 0.4, // Default: 0.4 m
              spacer_porosity = 0.79, // Default: 79%
              gap_spacer_porosity = 0.84, // Default: 84%
              air_gap_thickness = 0.8e-3, // Default: 0.8 mm
              wall_thickness = 62.0e-6, // Default: 62 microns
              polymer_conductivity = 0.35, // Default: 0.35 W/mK
              spacer_conductivity = 0.27, // Default: 0.27 W/mK, source: https://doi.org/10.1016/j.compositesa.2003.11.005
              wall_conductivity = 0.35; // Default: 0.35 W/mK
    PetscInt number_channels = 6; // Default: 6

    PetscOptionsGetReal(NULL, NULL, "-membrane_area", &membrane_area, NULL);
    PetscOptionsGetReal(NULL, NULL, "-membrane_thickness", &membrane_thickness, NULL);
    PetscOptionsGetReal(NULL, NULL, "-membrane_porosity", &membrane_porosity, NULL);
    PetscOptionsGetReal(NULL, NULL, "-pore_diameter", &pore_diameter, NULL);
    PetscOptionsGetReal(NULL, NULL, "-feed_channel_height", &feed_channel_height, NULL);
    PetscOptionsGetReal(NULL, NULL, "-cold_channel_height", &cool_channel_height, NULL);
    PetscOptionsGetReal(NULL, NULL, "-channel_width", &channel_width, NULL);
    PetscOptionsGetReal(NULL, NULL, "-spacer_porosity", &spacer_porosity, NULL);
    PetscOptionsGetReal(NULL, NULL, "-gap_spacer_porosity", &gap_spacer_porosity, NULL);
    PetscOptionsGetReal(NULL, NULL, "-air_gap_thickness", &air_gap_thickness, NULL);
    PetscOptionsGetReal(NULL, NULL, "-wall_thickness", &wall_thickness, NULL);
    PetscOptionsGetReal(NULL, NULL, "-polymer_conductivity", &polymer_conductivity, NULL);
    PetscOptionsGetReal(NULL, NULL, "-spacer_conductivity", &spacer_conductivity, NULL);
    PetscOptionsGetReal(NULL, NULL, "-wall_conductivity", &wall_conductivity, NULL);
    PetscOptionsGetInt(NULL, NULL, "-number_channels", &number_channels, NULL);

    dessal_data.membrane_area = membrane_area;
    dessal_data.membrane_thickness = membrane_thickness;
    dessal_data.membrane_porosity = membrane_porosity;
    dessal_data.pore_diameter = pore_diameter;
    dessal_data.feed_channel_height = feed_channel_height;
    dessal_data.cool_channel_height = cool_channel_height;
    dessal_data.channel_width = channel_width;
    dessal_data.number_channels = number_channels;
    dessal_data.spacer_porosity = spacer_porosity;
    dessal_data.gap_spacer_porosity = gap_spacer_porosity;
    dessal_data.air_gap_thickness = air_gap_thickness;
    dessal_data.wall_thickness = wall_thickness;
    dessal_data.polymer_conductivity = polymer_conductivity;
    dessal_data.spacer_conductivity = spacer_conductivity;
    dessal_data.wall_conductivity = wall_conductivity;

    // Iterative data
    PetscReal out_temperature_feed = entry_temperature_feed,
              out_temperature_cool = entry_temperature_cool,
              feed_membrane_temperature = entry_temperature_feed,
              gap_membrane_temperature = entry_temperature_feed,
              film_boundary_temperature = entry_temperature_cool,
              film_wall_temperature = entry_temperature_cool,
              cool_wall_temperature = entry_temperature_cool,
              out_salinity_feed = entry_salinity_feed,
              mass_flux = 0.0,
              heat_flux = 0.0,
              vapor_heat_flux = 0.0,
              feed_outflow_rate = feed_mass_flow_rate;

    dessal_data.out_temperature_feed = out_temperature_feed;
    dessal_data.out_temperature_cool = out_temperature_cool;
    dessal_data.feed_membrane_temperature = feed_membrane_temperature;
    dessal_data.gap_membrane_temperature = gap_membrane_temperature;
    dessal_data.film_boundary_temperature = film_boundary_temperature;
    dessal_data.film_wall_temperature = film_wall_temperature;
    dessal_data.cool_wall_temperature = cool_wall_temperature;
    dessal_data.out_salinity_feed = out_salinity_feed;
    dessal_data.mass_flux = mass_flux;
    dessal_data.heat_flux = heat_flux;
    dessal_data.vapor_heat_flux = vapor_heat_flux;
    dessal_data.feed_outflow_rate = feed_outflow_rate;

    // Aggregating all data
    entry_data->dessal_data = dessal_data;

    return 0;
}