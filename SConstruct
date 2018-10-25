import os
import sys
#import commands
import arch

vars = Variables()
vars.AddVariables(
  EnumVariable( 'order',
                'convergence order of the ADER-DG method',
                'none',
                allowed_values=('none',) + tuple(str(i) for i in range(2,33))
              ),
  PathVariable( 'buildDir', 'where to build the code', 'build', PathVariable.PathIsDirCreate ),
  EnumVariable( 'compileMode', 'mode of the compilation', 'release',
                allowed_values=('debug', 'release')
              ),
  EnumVariable( 'compiler',
                'Select the compiler (default: intel)',
                'intel',
                allowed_values=('intel', 'gcc')),
  EnumVariable( 'arch',
                'precision -- s for single- and d for double precision -- and architecture used. Warning: \'noarch\' calls the fall-back code and is outperformed by architecture-specific optimizations (if available) greatly.',
                'dnoarch',
                allowed_values=arch.getArchitectures()
              ),
)

# set environment
env = Environment(variables=vars)
env['ENV'] = os.environ

# generate help text
Help(vars.GenerateHelpText(env))

# handle unknown, maybe misspelled variables
unknownVariables = vars.UnknownVariables()

# exit in the case of unknown variables
if unknownVariables:
  ConfigurationError("*** The following build variables are unknown: " + str(unknownVariables.keys()))

if env['order'] == 'none':
  ConfigurationError("*** Convergence order not set.")

#
# preprocessor, compiler and linker
#

# Basic compiler setting
if env['compiler'] == 'intel':
    env['CC'] = 'icc'
    env['CXX'] = 'icpc'
    env['F90'] = 'ifort'
elif env['compiler'] == 'gcc':
    env['CC'] = 'gcc'
    env['CXX'] = 'g++'
    env['F90'] = 'gfortran'
else:
    assert(false)

#
# Common settings
#

# enforce restrictive C/C++-Code
env.Append(CFLAGS   = ['-Wall', '-Werror', '-ansi', '-qopenmp'],
           CXXFLAGS = ['-Wall', '-Werror', '-ansi', '-qopenmp'],
	   LINKFLAGS = ['-qopenmp'])

archFlags = arch.getFlags(env['arch'], env['compiler'])
env.Append( CFLAGS    = archFlags,
            CXXFLAGS  = archFlags,
            F90FLAGS  = archFlags,
            LINKFLAGS = archFlags )
env.Append(CPPDEFINES=['ALIGNMENT=' + str(arch.getAlignment(env['arch']))])

#
# Compile mode settings
#

# set (pre-)compiler flags for the compile modes
if env['compileMode'] == 'debug':
  env.Append(CFLAGS  = ['-O0', '-g'],
             CXXFLAGS = ['-O0', '-g'])
elif env['compileMode'] == 'release':
  env.Append(CPPDEFINES = ['NDEBUG'])
  env.Append(CFLAGS   = ['-O3'],
             CXXFLAGS = ['-O3'])

#
# Basic preprocessor defines
#

env.Append(CPPDEFINES=['CONVERGENCE_ORDER='+env['order']])
env.Append(CPPPATH=['#/src', '#/submodules/yateto/include'])
env.Append(CXXFLAGS=['--std=c++11'])

#
# setup the program name and the build directory
#
env['programName'] = 'lina_' + str(env['order'])
env['programFile'] = '%s/%s' %(env['buildDir'], env['programName'])

# build directory
env['buildDir'] = '%s/build_%s' %(env['buildDir'], env['programName'])

# get the source files
env.sourceFiles = []

Export('env')
SConscript('generated_code/SConscript', variant_dir='#/'+env['buildDir'], src_dir='#/', duplicate=0)
SConscript('src/SConscript', variant_dir='#/'+env['buildDir'], src_dir='#/', duplicate=0)
Import('env')

# build standard version
env.Program('#/'+env['programFile'], env.sourceFiles)
