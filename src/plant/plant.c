#include "solver.h"
#include "output.h"
#include "../dessal/dessal.h"

PetscErrorCode InitialGuess(Vec x, SolverCtx *solver_ctx)
{
    EntryData entry_data = solver_ctx->entry_data;
    DessalData dessal_data = entry_data.dessal_data;
    DM da = solver_ctx->da;
    PetscScalar *x_array;

    DMDAVecGetArray(da, x, &x_array);

    x_array[0] = dessal_data.out_temperature_feed;
    x_array[1] = dessal_data.out_temperature_cool;
    x_array[2] = dessal_data.feed_membrane_temperature;
    x_array[3] = dessal_data.gap_membrane_temperature;
    x_array[4] = dessal_data.film_boundary_temperature;
    x_array[5] = dessal_data.film_wall_temperature;
    x_array[6] = dessal_data.cool_wall_temperature;
    x_array[7] = dessal_data.out_salinity_feed;
    x_array[8] = dessal_data.mass_flux;
    x_array[9] = dessal_data.heat_flux;
    x_array[10] = dessal_data.vapor_heat_flux;
    x_array[11] = dessal_data.feed_outflow_rate;

    DMDAVecRestoreArray(da, x, &x_array);

    return 0;
}

PetscErrorCode PlantBalances(SNES snes, Vec x, Vec f, void *ctx)
{
    SolverCtx *solver_ctx = (SolverCtx *)ctx;
    EntryData entry_data = solver_ctx->entry_data;
    DessalData dessal_data = entry_data.dessal_data;
    DM da = solver_ctx->da;
    PetscScalar *x_array, *f_array;
    Vec x_local;

    DMGetLocalVector(da, &x_local);
    DMGlobalToLocal(da, x, INSERT_VALUES, x_local);
    DMDAVecGetArray(da, x_local, &x_array);
    DMDAVecGetArray(da, f, &f_array);

    //-----------------------------------------------------------------------------------------------------------------------------------------------//
    // Desalination module                                                                                                                           //
    //-----------------------------------------------------------------------------------------------------------------------------------------------//

    // Setting iterative data
    dessal_data.out_temperature_feed = x_array[0];
    dessal_data.out_temperature_cool = x_array[1];
    dessal_data.feed_membrane_temperature = x_array[2];
    dessal_data.gap_membrane_temperature = x_array[3];
    dessal_data.film_boundary_temperature = x_array[4];
    dessal_data.film_wall_temperature = x_array[5];
    dessal_data.cool_wall_temperature = x_array[6];
    dessal_data.out_salinity_feed = x_array[7];
    dessal_data.mass_flux = x_array[8];
    dessal_data.heat_flux = x_array[9];
    dessal_data.vapor_heat_flux = x_array[10];
    dessal_data.feed_outflow_rate = x_array[11];

    // Updating iterative data
    DessalBalance(&dessal_data);

    f_array[0] = x_array[0] - dessal_data.out_temperature_feed;
    f_array[1] = x_array[1] - dessal_data.out_temperature_cool;
    f_array[2] = x_array[2] - dessal_data.feed_membrane_temperature;
    f_array[3] = x_array[3] - dessal_data.gap_membrane_temperature;
    f_array[4] = x_array[4] - dessal_data.film_boundary_temperature;
    f_array[5] = x_array[5] - dessal_data.film_wall_temperature;
    f_array[6] = x_array[6] - dessal_data.cool_wall_temperature;
    f_array[7] = x_array[7] - dessal_data.out_salinity_feed;
    f_array[8] = x_array[8] - dessal_data.mass_flux;
    f_array[9] = x_array[9] - dessal_data.heat_flux;
    f_array[10] = x_array[10] - dessal_data.vapor_heat_flux;
    f_array[11] = x_array[11] - dessal_data.feed_outflow_rate;

    DMDAVecRestoreArray(da, x_local, &x_array);
    DMDAVecRestoreArray(da, f, &f_array);
    DMRestoreLocalVector(da, &x_local);

    return 0;
}

PetscErrorCode RunPlant(EntryData *entry_data)
{
    PetscFunctionBeginUser;

    SolverCtx solver_ctx;
    Vec solution;
    Mat jac;
    char file[256] = "./results/report.csv";

    SolverCtxBuild(&solver_ctx, entry_data);

    DMCreateGlobalVector(solver_ctx.da, &solution);
    DMCreateMatrix(solver_ctx.da, &jac);

    InitialGuess(solution, &solver_ctx);

    SNESSetFunction(solver_ctx.snes, NULL, PlantBalances, &solver_ctx);
    SNESSetJacobian(solver_ctx.snes, jac, jac, SNESComputeJacobianDefault, NULL);

    SNESSolve(solver_ctx.snes, NULL, solution);

    ExportToFile(&solution, entry_data, file);

    VecDestroy(&solution);
    MatDestroy(&jac);
    SolverCtxDestroy(&solver_ctx);

    return 0;
}