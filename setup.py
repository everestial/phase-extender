import setuptools
setuptools.setup(
    name='phase-Extender',
    version='1.0',
    scripts=['./scripts/myscript'],
    author='Me',
    description='This runs my script which is great.',
    packages=['lib.myscript'],
    install_requires=[
        'setuptools',
        'pandas >= 0.22.0',
        'numpy >= 1.16.0'
    ],
    python_requires='>=3.5'
)