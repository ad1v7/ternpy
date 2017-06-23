from setuptools import setup, find_packages


def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='ternpy',
      version='0.2.0',
      long_description=readme(),
      description='Ternary Plot Tool',
      url='http://github.com/ad1v7/ternpy',
      author=['Marcin Kirsz', 'Sebastiaan van de Bund'],
      author_email=['marcin.kirsz@gmail.com', 'sebvdb@me.com'],
      license='MIT',
      packages=find_packages(),
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'sympy',
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      include_package_data=True,
      entry_points={
          'console_scripts': ['ternpy=ternpy.cli:main'],
      },
      zip_safe=False)
