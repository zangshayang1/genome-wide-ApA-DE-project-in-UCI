from distutils.core import setup

setup(
      name = 'APAP_master',
      version = '0.1dev',
      description = 'Please refer README.txt',
      author = 'Shayang Zang',
      author_email = 'szang@uci.edu',
      url = 'https://github.com/zangshayang1/genome-wide-ApA-DE-project-in-UCI',
      license = open('LICENCE.txt').read(),
      packages = ['APAP',],
      long_description=open('README.txt').read(),
      py_modules = ['io',
                    'const',
                    'logger',
                    'reference',
                    'utils',
                    'myclasses',
                    'filter',
                    'cleaner',
                    'counter',
                    'merger',
                    'summer',
                    'rc2pcter',
                    'top2selecter',
                    'redcalculator',
		            'hister'],
      scripts = ['apap.py']
        )
