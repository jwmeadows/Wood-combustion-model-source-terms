/* This file generated automatically. */ 
/* Do not modify. */ 
#include "udf.h" 
#include "prop.h" 
#include "dpm.h" 
extern DEFINE_INIT(density_initial, domain);
extern DEFINE_PROPERTY(rho_wood, c, t);   /* actual wood density profile computation*/
extern DEFINE_PROPERTY(thermal_conductivity, c, t);
extern DEFINE_EXECUTE_AT_END(execute_at_end); /*print density difference at the end of each time step */
extern DEFINE_EXECUTE_AT_END(execute_at_end_solid); /* loop over solid to assign char distribution values */
extern DEFINE_SOURCE(oxygen_sink, c, t, ds, eqn);
extern DEFINE_SOURCE(co2_flux, c, t, ds, eqn);
extern DEFINE_SOURCE(co_flux, c, t, ds, eqn);
extern DEFINE_SOURCE(energy_source, c, t, ds, eqn);
extern DEFINE_EXECUTE_ON_LOADING(load, libname);
extern DEFINE_SOURCE(CO_kinetics, c, t, ds, eqn);
extern DEFINE_SOURCE(CO2_kinetics, c, t, ds, eqn);
extern DEFINE_SOURCE(CH4_kinetics, c, t, ds, eqn);
extern DEFINE_SOURCE(hemicellulose_hod, c, t, ds, eqn);
extern DEFINE_SOURCE(cellulose_hod, c, t, ds, eqn);
extern DEFINE_SOURCE(lignin_hod, c, t, ds, eqn);
__declspec(dllexport) UDF_Data udf_data[] = { 
{"density_initial", (void (*)(void))density_initial, UDF_TYPE_INIT},
{"rho_wood", (void (*)(void))rho_wood, UDF_TYPE_PROPERTY},   /* actual wood density profile computation*/
{"thermal_conductivity", (void (*)(void))thermal_conductivity, UDF_TYPE_PROPERTY},
{"execute_at_end", (void (*)(void))execute_at_end, UDF_TYPE_EXECUTE_AT_END}, /*print density difference at the end of each time step */
{"execute_at_end_solid", (void (*)(void))execute_at_end_solid, UDF_TYPE_EXECUTE_AT_END}, /* loop over solid to assign char distribution values */
{"oxygen_sink", (void (*)(void))oxygen_sink, UDF_TYPE_SOURCE},
{"co2_flux", (void (*)(void))co2_flux, UDF_TYPE_SOURCE},
{"co_flux", (void (*)(void))co_flux, UDF_TYPE_SOURCE},
{"energy_source", (void (*)(void))energy_source, UDF_TYPE_SOURCE},
{"load", (void (*)(void))load, UDF_TYPE_EXECUTE_ON_LOADING},
{"CO_kinetics", (void (*)(void))CO_kinetics, UDF_TYPE_SOURCE},
{"CO2_kinetics", (void (*)(void))CO2_kinetics, UDF_TYPE_SOURCE},
{"CH4_kinetics", (void (*)(void))CH4_kinetics, UDF_TYPE_SOURCE},
{"hemicellulose_hod", (void (*)(void))hemicellulose_hod, UDF_TYPE_SOURCE},
{"cellulose_hod", (void (*)(void))cellulose_hod, UDF_TYPE_SOURCE},
{"lignin_hod", (void (*)(void))lignin_hod, UDF_TYPE_SOURCE},
}; 
__declspec(dllexport) int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data); 
#include "version.h" 
__declspec(dllexport) void UDF_Inquire_Release(int *major, int *minor, int *revision) 
{ 
*major = RampantReleaseMajor; 
*minor = RampantReleaseMinor; 
*revision = RampantReleaseRevision; 
} 
