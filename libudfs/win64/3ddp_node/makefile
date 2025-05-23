###########################################################
##
## Copyright 1987-2023 ANSYS, Inc. All Rights Reserved. 
## All Rights Reserved  
##
###########################################################
# 
# purpose:
#    This file is the makefile for NT UDF builds.
#
# usage:
#    Change file name 'makefile_nt.udf' to 'makefile'.
#    Place the makefile file in the target dir.
#    Place the helper file 'user_nt.udf' in same directory.
#    Edit user_nt.udf file.
#    Execute an 'nmake'  (no arguments needed).
#
# notes:
#    The user should not need to change this file.
#    The following macros are expected to be in the user_nt.udf file:
#      SOURCES, VERSION, PARALLEL_NODE and USER_OBJECTS.
#    The macro FLUENT_LIB is automatically defined here based on VERSION
#      and PARALLEL_NODE.
#    Only this makefile and the helper file 'user_nt.udf' needs to be 
#      the target dir, all other files are automatically generated, 
#      and removed in an 'nmake clean'.
#    The source file is expected in ..\..\src directory.
#    Make sure the VC compiler is set up for command line mode, use
#      the file Vcvars32.bat from Visual Studio to assist this. 
#    This makefile uses the public domain program 'sed', make sure
#      sed is in your executable path.
#    The '^' char in the SRC literal is a quote for the '\' char,
#      otherwise interpretted as a line continuation char.
#    Use the 'nmake clean' command to restart builds.
###########################################################

!INCLUDE user_nt.udf

# These need to be updated for each release.
LIB_RELEASE=2320
RELEASE=23.2.0

# If a local path can be found, use it
!IF (("$(FLUENT_INC_LOCAL_RESERVED)" != "") && ("$(FLUENT_INC)" == ""))
FLUENT_INC = $(FLUENT_INC_LOCAL_RESERVED)
!ENDIF

# it is assumed $(FLUENT_INC)is already set.
!IF !DEFINED(FLUENT_INC)
!  ERROR You need to define the environmental variable FLUENT_INC.
!ENDIF

!IF !DEFINED(FLUENT_ARCH)

!IF ("$(PROCESSOR_ARCHITECTURE)" == "ALPHA")
FLUENT_ARCH = ntalpha
!ENDIF

!IF ("$(PROCESSOR_ARCHITECTURE)" == "x86")
FLUENT_ARCH = ntx86
!ENDIF

!IF ("$(PROCESSOR_ARCHITECTURE)" == "AMD64")
FLUENT_ARCH = win64 
!ENDIF

!IF ("$(PROCESSOR_ARCHITEW6432)" == "AMD64")
!  MESSAGE You seem to be working on a 64-bit machine, but in a 32-bit build environment.
!  MESSAGE If you really want to build the UDF for the "ntx86" FLUENT_ARCH, then please
!  MESSAGE set the "FLUENT_ARCH" environment variable to "ntx86" before calling this.
!  ERROR   Aborting.
!ENDIF

!ENDIF

# set FLUENT_LIB based on VERSION macro

!IF (("$(VERSION)" == "2d") || ("$(VERSION)" == "3d") || \
     ("$(VERSION)" == "2ddp") || ("$(VERSION)" == "3ddp"))
FLUENT_LIB = fl$(LIB_RELEASE).lib
LIBS = /Libpath:"$(FLUENT_INC)"\fluent$(RELEASE)\$(FLUENT_ARCH)\$(VERSION)
!ENDIF

!IF (("$(VERSION)" == "2d_host") || ("$(VERSION)" == "3d_host") || \
     ("$(VERSION)" == "2ddp_host") || ("$(VERSION)" == "3ddp_host"))
COMM = net
FLUENT_LIB =  mport.lib fl$(LIB_RELEASE).lib
LIBS = /Libpath:"$(FLUENT_INC)"\fluent$(RELEASE)\$(FLUENT_ARCH)\$(VERSION)  /Libpath:"$(FLUENT_INC)"\fluent$(RELEASE)\multiport\$(FLUENT_ARCH)\$(COMM)\shared
!ENDIF

!IF (("$(VERSION)" == "2d_node") || ("$(VERSION)" == "3d_node") || \
     ("$(VERSION)" == "2ddp_node") || ("$(VERSION)" == "3ddp_node"))
!IF (("$(PARALLEL_NODE)" != "none") && ("$(PARALLEL_NODE)" != "net"))
PARALLEL_NODE = mpi
!ENDIF
COMM = mpi
FLUENT_LIB =  mport.lib fl_$(PARALLEL_NODE)$(LIB_RELEASE).lib
LIBS = /Libpath:"$(FLUENT_INC)"\fluent$(RELEASE)\$(FLUENT_ARCH)\$(VERSION)  /Libpath:"$(FLUENT_INC)"\fluent$(RELEASE)\multiport\$(FLUENT_ARCH)\$(COMM)\shared
!ENDIF


###########################################################

SRC = ..\..\src^\

CC = cl -nologo

# compile flags 
#  /c for compile only,  /Za for ansi c and __STDC__ is set
CFLAGS = /c /DUDF_EXPORTING /DUDF_NT /DWIN64 /EHa /wd4224

TARGET = libudf.dll
UDFDATA = udf_names.c

# generate object names from source names
SOURCES2 = $(CSOURCES:..\..\src\=)
SRC_OBJECT = $(SOURCES2:.cpp=.obj)
SRC_OBJECT = $(SRC_OBJECT:.c=.obj)
UDF_OBJECT = $(UDFDATA:.c=.obj) 

OBJECTS = $(UDF_OBJECT) $(SRC_OBJECT) $(USER_OBJECTS)

INCLUDES= -I. \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\$(FLUENT_ARCH)\$(VERSION) \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\main \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\addon-wrapper \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\io \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\species \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\pbns \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\numerics \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\sphysics \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\storage \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\mphase \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\bc \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\models \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\material \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\amg \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\util \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\mesh \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\udf \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\ht \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\dx \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\turbulence \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\acoustics \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\parallel \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\etc \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\ue \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\dpm \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\src\dbns \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\cortex\src \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\client\src \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\tgrid\src \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\multiport\src \
          -I"$(FLUENT_INC)"\fluent$(RELEASE)\multiport\mpi_wrapper\src \
          -I"$(FLUENT_INC)"\include

GPU_LIB=
!IF ("$(GPU_SUPPORT)" == "on")
INCLUDES = $(INCLUDES) \
           -I"$(FLUENT_INC)"\fluent$(RELEASE)\multiport\gpu_wrapper\include
GPU_LIB = OpenCL.lib
LIBS = $(LIBS) /Libpath:"$(FLUENT_INC)"\fluent$(RELEASE)\multiport\gpu_wrapper\$(FLUENT_ARCH)\stub
!ENDIF

default:   $(TARGET)

$(UDF_OBJECT): $(UDFDATA)
	$(CC) $(CFLAGS) $(INCLUDES)  $**
	
$(SRC_OBJECT): $(CSOURCES)
	@echo # Generating ud_io1.h
	resolve.exe -udf $(CSOURCES) -head_file ud_io1.h
	$(CC) $(CFLAGS) $(HOST_FLAGS) $(INCLUDES)  $**
	
$(TARGET): makefile user_nt.udf $(UDF_OBJECT) $(SRC_OBJECT) 
	@echo # Linking $@ because of $?
	link -nologo $(LIBS) /dll   /out:$(TARGET) \
	      $(OBJECTS)  $(FLUENT_LIB) $(GPU_LIB)
	
$(UDFDATA): makefile $(SRC_OBJECT)
	@echo # Generating $@ because of $?
	@echo /* This file generated automatically. */ > $@
	@echo /*          Do not modify.            */ >> $@
	@echo #include "udf.h" >> $@
	@echo #include "prop.h" >> $@
	@echo #include "dpm.h" >> $@
		
	@sed -n "s/^ *\(DEFINE_[_A-Z]*(.*)\)/extern \1;/p" $(CSOURCES) >> $@
	
	@echo __declspec(dllexport) UDF_Data udf_data[] = { >> $@
	
	@sed -n "s/^ *DEFINE_\([_A-Z]*\)( *\([_a-zA-Z0-9]*\).*)/\{\"\2\", (void (*)(void))\2, UDF_TYPE_\1\},/p" $(CSOURCES) >> $@
	
	@echo }; >> $@
	
	@echo __declspec(dllexport) int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data); >> $@
	@echo #include "version.h" >> $@
	@echo __declspec(dllexport) void UDF_Inquire_Release(int *major, int *minor, int *revision) >> $@
	@echo { >> $@
	@echo   *major = RampantReleaseMajor; >> $@
	@echo   *minor = RampantReleaseMinor; >> $@
	@echo   *revision = RampantReleaseRevision; >> $@
	@echo } >> $@
	

clean:
	@del  *.dll *.obj *.ilk *.exp *.lib *.pdb *.c 

tellfluentarch:
	@echo FLUENT_ARCH: $(FLUENT_ARCH)
	@echo PROCESSOR_ARCHITECTURE: $(PROCESSOR_ARCHITECTURE)
	@echo PROCESSOR_ARCHITEW6432: $(PROCESSOR_ARCHITEW6432)


$(OBJECTS): $(HSOURCES)

