from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

DESCRIPTION  = """
A set of scripts for working with Solexa data. Not (yet) intended for general distribution.
For any question, please ask Simon van Heeringen
"""

setup (name = 'solexatools',
		version = '0.22',
		description = DESCRIPTION,
		author='Simon van Heeringen',
		author_email='s.vanheeringen@ncmls.ru.nl',
		license='MIT',
		packages=['solexatools'],
		cmdclass = {'build_ext': build_ext},
		ext_modules = [Extension("solexatools.track", ["solexatools/track.pyx"])],
		scripts=[
			'scripts/assign_peaks.py', 
			'scripts/peakstats.py', 
			'scripts/peak_heap.py', 
			'scripts/solexa_analysis.py',
			'scripts/coverage_profile.py',
			'scripts/other_profile.py',
			'scripts/polII_index.py',
			'scripts/rpkm_rnapii.py',
			],
)
