#ifndef HYPRE_PARCSR_INTERFACE_INTERPRETER
#define HYPRE_PARCSR_INTERFACE_INTERPRETER

#include "interpreter.h"
#include "HYPRE_MatvecFunctions.h"

#ifdef __cplusplus
extern "C" {
#endif

int
HYPRE_StructSetupInterpreter( mv_InterfaceInterpreter *i );

int
HYPRE_StructSetupMatvec(HYPRE_MatvecFunctions * mv);


#ifdef __cplusplus
}
#endif

#endif
