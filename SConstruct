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
  BoolVariable( 'unitTests', 'Build unit tests', False),
  BoolVariable( 'libxsmm', 'Generate kernels with LIBXSMM.', True),
  BoolVariable( 'pspamm', 'Generate kernels with PSpaMM.', True),
  BoolVariable( 'sparse', 'Use sparse memory layout.', True)
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

# Link to MKL
env.Append(LIBS=['mkl_intel_lp64', 'mkl_sequential', 'mkl_core'])
mklroot = env['ENV']['MKLROOT']
if mklroot:
  env.Append(LIBPATH=os.path.join(mklroot, 'lib', 'intel64'))
  env.Append(CPPPATH=os.path.join(mklroot, 'include'))

# enforce restrictive C/C++-Code
env.Append(CFLAGS   = ['-Wall', '-ansi'],
           CXXFLAGS = ['-Wall', '-ansi'])

if env['compiler'] == 'intel':
  env.Append( CFLAGS   = ['-qopenmp'],
              CXXFLAGS = ['-qopenmp'],
              LINKFLAGS = ['-qopenmp'])
else:
  env.Append( CFLAGS   = ['-fopenmp'],
              CXXFLAGS = ['-fopenmp'],
              LINKFLAGS = ['-fopenmp'])

archFlags = arch.getFlags(env['arch'], env['compiler'])
env.Append( CFLAGS    = archFlags,
            CXXFLAGS  = archFlags,
            F90FLAGS  = archFlags,
            LINKFLAGS = archFlags )
env.Append(CPPDEFINES=arch.getDefines(env['arch']))

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
generators = list()
if env['libxsmm']:
  generators.append('x')
if env['pspamm']:
  generators.append('p')
ml = 'sp' if env['sparse'] else 'dn'
programSuffix = '_{}_{}_{}_{}'.format(env['arch'], ''.join(generators), ml, env['order'])
env['programName'] = 'lina' + programSuffix
env['programFile'] = '%s/%s' %(env['buildDir'], env['programName'])
unitTestProgramFile = os.path.join(env['buildDir'], 'unit_tests' + programSuffix)

# build directory
env['buildDir'] = '%s/build_%s' %(env['buildDir'], env['programName'])

# get the source files
env.sourceFiles = []
env.generatedTestSourceFiles = []

Export('env')
SConscript('generated_code/SConscript', variant_dir='#/'+env['buildDir'], src_dir='#/', duplicate=0)
SConscript('src/SConscript', variant_dir='#/'+env['buildDir'], src_dir='#/', duplicate=0)
Import('env')

# build standard version
env.Program('#/'+env['programFile'], env.sourceFiles)

if env['unitTests']: 
  env.Tool('cxxtest')
  sourceFiles = filter(lambda sf: os.path.basename(str(sf[0])) != 'main.o', env.sourceFiles)
  env.CxxTest(unitTestProgramFile, list(sourceFiles) + env.generatedTestSourceFiles)
