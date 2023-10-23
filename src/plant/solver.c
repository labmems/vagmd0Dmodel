#include "solver.h"

#define NUM_VAR 12

PetscErrorCode SolverCtxBuild(SolverCtx *solver_ctx, EntryData *entry_data)
{
    PetscFunctionBeginUser;
    SNES snes;
    SNESLineSearch snesls;
    KSP ksp;
    DM da;

    SNESCreate(PETSC_COMM_WORLD, &snes);
    SNESSetType(snes, SNESNEWTONLS);
    SNESGetLineSearch(snes, &snesls);
    SNESLineSearchSetType(snesls, SNESLINESEARCHL2);
    SNESLineSearchSetTolerances(snesls, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1000);
    SNESSetTolerances(snes, 1.0e-10, 1.0e-10, PETSC_DEFAULT, 2000, -1);
    SNESGetKSP(snes, &ksp);
    KSPSetTolerances(ksp, 1.0e-12, 1.0e-12, PETSC_DEFAULT, 2000);
    KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
    SNESSetFromOptions(snes);

    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, NUM_VAR, 1, 1, NULL, &da);
    DMSetUp(da);

    solver_ctx->snes = snes;
    solver_ctx->da = da;
    solver_ctx->entry_data = *entry_data;

    return 0;
}

PetscErrorCode SolverCtxDestroy(SolverCtx *solver_ctx)
{
    PetscFunctionBeginUser;
    SNESDestroy(&solver_ctx->snes);
    DMDestroy(&solver_ctx->da);

    return 0;
}