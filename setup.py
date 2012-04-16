from distutils.core import setup, Extension

DESCRIPTION  = """
p53scan (http://www.ncmls.nl/bioinfo/p53scan/) is an algorithm to 
locate p53 binding sites in DNA sequences. It is described in: 
Smeenk L, van Heeringen SJ, Koeppel M, Driel MA, Bartels SJ, 
Akkers RC, Denissov S, Stunnenberg HG, Lohrum M. Characterization
of genome-wide p53-binding sites upon stress response. 
Nucleic Acids Res. (2008)"""

module1 = Extension('PftScan.core',
					sources = ['PftScan/coremodule.c'])

setup (name = 'p53scan',
		version = '1.05',
		description = DESCRIPTION,
		author='Simon van Heeringen',
		author_email='s.vanheeringen@ncmls.ru.nl',
		url='http://www.ncmls.nl/bioinfo/p53scan',
		license='MIT',
		ext_modules = [module1],
		packages=['PftScan'],
		scripts=['scripts/p53scan.py', 'scripts/pwmscan.py', 'scripts/scantriplet.py', 'scripts/p53scan_gui.py', 'scripts/windows_postinstall.py', 'scripts/p63scan.py', 'scripts/p63scan_gui.py'],
)
