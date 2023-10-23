#include "properties.h"
#include "../entrydata/entrydata.h"

/*
Salt water thermophysical propeties

Reference: K.G. Nayar, M.H. Sharqawy, L.D. Banchik, J.H. Lienhard IV, Thermophysical properties of seawater: A review and new correlations that
           include pressure dependence. Desalination 390 (2016) 1-24. https://doi.org/10.1016/j.desal.2016.02.024
*/

PetscReal SaltWaterDensity(PetscReal temperature, PetscReal salinity)
{
    PetscReal a[5] = {9.999e2,
                      2.034e-2,
                      -6.162e-3,
                      2.261e-5,
                      -4.657e-8};
    PetscReal b[5] = {8.020e2,
                      -2.001,
                      1.677e-2,
                      -3.060e-5,
                      -1.613e-5};
    PetscReal temperature_part, salinity_part;

    temperature_part = a[0] + a[1] * temperature;
    temperature_part += a[2] * temperature * temperature;
    temperature_part += a[3] * temperature * temperature * temperature;
    temperature_part += a[4] * temperature * temperature * temperature * temperature;

    salinity_part = b[0] * salinity;
    salinity_part += b[1] * salinity * temperature;
    salinity_part += b[2] * salinity * temperature * temperature;
    salinity_part += b[3] * salinity * temperature * temperature * temperature;
    salinity_part += b[4] * salinity * salinity * temperature * temperature;

    return temperature_part + salinity_part;
}

PetscReal SaltWaterSpecificHeat(PetscReal temperature, PetscReal salinity)
{
    PetscReal a[3] = {5328.0,
                      -9.76e1,
                      4.04e-1};
    PetscReal b[3] = {-6.913,
                      7.351e-1,
                      -3.15e-3};
    PetscReal c[3] = {9.6e-3,
                      -1.927e-3,
                      8.23e-6};
    PetscReal d[3] = {2.5e-6,
                      1.666e-6,
                      -7.125e-9};
    PetscReal alt_salinity, alt_salinity2, abs_temperature, A, B, C, D, specific_heat;

    alt_salinity = 1000.0 * salinity;
    alt_salinity2 = alt_salinity * alt_salinity;

    A = a[0] + a[1] * alt_salinity + a[2] * alt_salinity2;
    B = b[0] + b[1] * alt_salinity + b[2] * alt_salinity2;
    C = c[0] + c[1] * alt_salinity + c[2] * alt_salinity2;
    D = d[0] + d[1] * alt_salinity + d[2] * alt_salinity2;

    abs_temperature = temperature + 273.15;

    specific_heat = A + B * abs_temperature + C * abs_temperature * abs_temperature;
    specific_heat += D * abs_temperature * abs_temperature * abs_temperature;

    return specific_heat;
}

// Exceptionally taken from https://doi.org/10.5004/dwt.2010.1079
PetscReal SaltWaterDynViscosity(PetscReal temperature, PetscReal salinity)
{
    PetscReal a[3] = {0.0428,
                      0.00123,
                      0.000131};
    PetscReal b[3] = {-0.03724,
                      0.01859,
                      -0.00271};
    PetscReal c[3] = {4.2844e-5,
                      0.157,
                      -91.296};
    PetscReal pure_viscosity, viscosity, alt_temperature, alt_salinity, ionic_strength, ionic_strength2, ionic_strength3;

    alt_salinity = salinity / 1.00472;
    alt_temperature = temperature + 64.993;
    ionic_strength = 19.915 * alt_salinity / (1.0 - 1.00487 * alt_salinity);
    ionic_strength2 = ionic_strength * ionic_strength;
    ionic_strength3 = ionic_strength * ionic_strength2;

    pure_viscosity = c[1] * alt_temperature * alt_temperature + c[2];
    pure_viscosity = c[0] + 1.0 / pure_viscosity;

    viscosity = b[0] * ionic_strength + b[1] * ionic_strength2 + b[2] * ionic_strength3;
    viscosity *= PetscLog10Real(1000.0 * pure_viscosity);
    viscosity += a[0] * ionic_strength + a[1] * ionic_strength2 + a[2] * ionic_strength3;
    viscosity = pure_viscosity * PetscPowReal(10.0, viscosity);

    return viscosity;
}

PetscReal SaltWaterThermalConductivity(PetscReal temperature, PetscReal salinity)
{
    PetscReal b[4] = {0.797015,
                      -0.251242,
                      0.096437,
                      -0.032696};
    PetscReal alt_salinity, dimless_temperature, thermal_conductivity;

    alt_salinity = 1000.0 * salinity;
    dimless_temperature = (temperature + 273.15) / 300.0;

    thermal_conductivity = b[0] * PetscPowReal(dimless_temperature, -0.194);
    thermal_conductivity += b[1] * PetscPowReal(dimless_temperature, -4.717);
    thermal_conductivity += b[2] * PetscPowReal(dimless_temperature, -6.385);
    thermal_conductivity += b[3] * PetscPowReal(dimless_temperature, -2.134);
    thermal_conductivity /= (1.0 + 0.00022 * alt_salinity);

    return thermal_conductivity;
}

// The vapor pressure framework was conducted using https://doi.org/10.1016/j.ijheatmasstransfer.2013.07.051
PetscReal PureVaporPressure(PetscReal temperature)
{
    PetscReal vapor_pressure;

    vapor_pressure = 23.1964 - 3816.44 / (temperature + 227.02);
    vapor_pressure = PetscExpReal(vapor_pressure);

    return vapor_pressure;
}

PetscReal ActivityCoefficient(PetscReal salinity)
{
    PetscReal activity_coefficient, molar_fraction;

    molar_fraction = water_molar_mass * salinity / ((1.0 - salinity) * salt_molar_mass + salinity * water_molar_mass);

    activity_coefficient = 1.0 - 0.5 * molar_fraction - 10.0 * molar_fraction * molar_fraction;
    activity_coefficient *= (1.0 - molar_fraction);

    return activity_coefficient;
}

PetscReal VaporPressure(PetscReal temperature, PetscReal salinity)
{
    PetscReal pure_vapor_pressure, activity_coefficient;

    pure_vapor_pressure = PureVaporPressure(temperature);
    activity_coefficient = ActivityCoefficient(salinity);

    return pure_vapor_pressure * activity_coefficient;
}

// Exceptionally taken from https://doi.org/10.5004/dwt.2010.1079
PetscReal SaltWaterLatentHeat(PetscReal temperature, PetscReal salinity)
{
    PetscReal a[5] = {2.501e6,
                      -2.369e3,
                      2.678e-1,
                      -8.103e-3,
                      -2.079e-5};
    PetscReal pure_latent_heat, latent_heat_vaporization;

    pure_latent_heat = a[0] + a[1] * temperature;
    pure_latent_heat += a[2] * temperature * temperature;
    pure_latent_heat += a[3] * temperature * temperature * temperature;
    pure_latent_heat += a[4] * temperature * temperature * temperature * temperature;

    latent_heat_vaporization = pure_latent_heat * (1.0 - salinity);

    return latent_heat_vaporization;
}

PetscErrorCode SaltWaterPropBuild(SaltWaterProperties *salt_water_prop, PetscReal temperature, PetscReal salinity)
{
    PetscFunctionBeginUser;

    PetscReal density, specific_heat, dyn_viscosity, thermal_conductivity, vapor_pressure, latent_heat_vaporization;

    density = SaltWaterDensity(temperature, salinity);
    specific_heat = SaltWaterSpecificHeat(temperature, salinity);
    dyn_viscosity = SaltWaterDynViscosity(temperature, salinity);
    thermal_conductivity = SaltWaterThermalConductivity(temperature, salinity);
    vapor_pressure = VaporPressure(temperature, salinity);
    latent_heat_vaporization = SaltWaterLatentHeat(temperature, salinity);

    salt_water_prop->density = density;
    salt_water_prop->specific_heat = specific_heat;
    salt_water_prop->dyn_viscosity = dyn_viscosity;
    salt_water_prop->thermal_conductivity = thermal_conductivity;
    salt_water_prop->prandtl = dyn_viscosity * specific_heat / thermal_conductivity;
    salt_water_prop->vapor_pressure = vapor_pressure;
    salt_water_prop->latent_heat_vaporization = latent_heat_vaporization;

    return 0;
}

/*
Moist air properties

Reference: P.T. Tsilingiris, Review and critical comparative evaluation of moist air thermophysical properties at the temperature range between
           0 and 100 °C for Engineering Calculations. Renew. and Sust. Energy Reviews 83 (2018) 50-63. https://doi.org/10.1016/j.rser.2017.10.072

Assumption: Air saturated with water vapor (RH = 100%)
*/

PetscReal MoistAirDensity(PetscReal temperature)
{
    PetscReal sd[4] = {1.293393662,
                       -5.538444326e-3,
                       3.860201577e-5,
                       -5.2536065e-7};
    PetscReal density;

    density = sd[0] + sd[1] * temperature;
    density += sd[2] * temperature * temperature;
    density += sd[3] * temperature * temperature * temperature;

    return density;
}

PetscReal MoistAirSpecificHeat(PetscReal temperature)
{
    PetscReal sc[6] = {1.00457142,
                       2.05063275e-3,
                       -1.6315370e-4,
                       6.2123003e-6,
                       -8.8304788e-8,
                       5.07130703e-10};
    PetscReal specific_heat;

    specific_heat = sc[0] + sc[1] * temperature;
    specific_heat += sc[2] * temperature * temperature;
    specific_heat += sc[3] * temperature * temperature * temperature;
    specific_heat += sc[4] * temperature * temperature * temperature * temperature;
    specific_heat += sc[5] * temperature * temperature * temperature * temperature * temperature;
    specific_heat *= 1000.0;

    return specific_heat;
}

PetscReal MoistAirDynViscosity(PetscReal temperature)
{
    PetscReal sv[5] = {1.715747771e-5,
                       4.722402075e-8,
                       -3.663027156e-10,
                       1.873236686e-12,
                       -8.050218737e-14};
    PetscReal dyn_viscosity;

    dyn_viscosity = sv[0] + sv[1] * temperature;
    dyn_viscosity += sv[2] * temperature * temperature;
    dyn_viscosity += sv[3] * temperature * temperature * temperature;
    dyn_viscosity += sv[4] * temperature * temperature * temperature * temperature;

    return dyn_viscosity;
}

PetscReal MoistAirThermalConductivity(PetscReal temperature)
{
    PetscReal sk[5] = {24.0073953e-3,
                       7.278410162e-5,
                       -1.788037411e-7,
                       -1.351703529e-9,
                       -3.322412767e-11};
    PetscReal thermal_conductivity;

    thermal_conductivity = sk[0] + sk[1] * temperature;
    thermal_conductivity += sk[2] * temperature * temperature;
    thermal_conductivity += sk[3] * temperature * temperature * temperature;
    thermal_conductivity += sk[4] * temperature * temperature * temperature * temperature;

    return thermal_conductivity;
}

PetscErrorCode MoistAirPropBuild(MoistAirProperties *moist_air_prop, PetscReal temperature)
{
    PetscReal density, specific_heat, dyn_viscosity, thermal_conductivity;

    density = MoistAirDensity(temperature);
    specific_heat = MoistAirSpecificHeat(temperature);
    dyn_viscosity = MoistAirDynViscosity(temperature);
    thermal_conductivity = MoistAirThermalConductivity(temperature);

    moist_air_prop->density = density;
    moist_air_prop->specific_heat = specific_heat;
    moist_air_prop->dyn_viscosity = dyn_viscosity;
    moist_air_prop->thermal_conductivity = thermal_conductivity;
    moist_air_prop->prandtl = dyn_viscosity * specific_heat / thermal_conductivity;

    return 0;
}