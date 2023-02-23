from setuptools import setup, find_packages

MAJOR = 4
MINOR = 10
PATCH = 2
VERSION = "{}.{}.{}".format(MAJOR, MINOR, PATCH)

with open("ttclust/version.py", "w") as f:
    f.write("__version__ = '{}'\n".format(VERSION))


setup(
    name='ttclust',
    version=VERSION,
    url='https://github.com/tubiana/TTClust',
    license='GPL3',
    author='Thibault Tubiana',
    author_email='tubiana.thibault@gmail.com',
    description='TTclust : A molecular simulation clustering program',
    platforms=["Linux", "Solaris", "Mac OS-X", "darwin", "Unix", "win32"],
    setup_requires = ['cython'],
    install_requires=['argparse',
                      'argcomplete',
                      'cython',
                      'progressbar2',
                      'matplotlib',
                      'numpy',
                      'prettytable',
                      'pandas',
                      'scipy >= 0.18',
                      'scikit-learn',
                      'numba',
                      'hashlib'
                      'mdtraj >= 1.7'],

    entry_points={'console_scripts':['ttclust=ttclust.ttclust:main']},
        #'gui_scripts':['ttclustGUI=ttclust.ttclustGUI:main'],





    packages=find_packages(),
)
