from setuptools import setup, find_packages

setup(
    name='ttclust',
    version='4.6',
    url='https://github.com/tubiana/TTClust',
    license='GPL3',
    author='Thibault Tubiana',
    author_email='tubiana.thibault@gmail.com',
    description='A molecular simulation clustering program',
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
                      'sklearn ',
                      'mdtraj >= 1.7',
                      'msgpack',
                      'RXPY>=0.1.0',
                      'wxpython==4.0.0b1',
                      'Pillow==4.3.0',
                      'psutil==5.4.2',
                      'gooey'],
    entry_points={
        'console_scripts':['ttclust=ttclust.ttclust:main'],
        'gui_scripts':['ttclustGUI=ttclust.ttclustGUI:main'],
    },





    packages=find_packages(),
)
