from distutils.core import setup

setup(
      name = 'APAP_master',
      version = '0.1dev',
      description = 'Placeholder',
      author = 'Shayang Zang',
      author_email = 'szang@uci.edu',
      url = 'Placeholder',
      license = open('LICENCE.rtf').read(),
      packages = ['APAP',],
      long_description=open('README.rtf').read(),
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
