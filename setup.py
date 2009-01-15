from distutils.core import setup

DESCRIPTION  = """
A set of scripts for working with Solexa data. Not (yet) intended for general distribution.
For any question, please ask Simon van Heeringen
"""

setup (name = 'SolexaTools',
		version = '0.1a',
		description = DESCRIPTION,
		author='Simon van Heeringen',
		author_email='s.vanheeringen@ncmls.ru.nl',
		license='MIT',
		packages=['SolexaTools'],
		scripts=['scripts/assign_peaks.py', 'scripts/peakstats.py'],
)
