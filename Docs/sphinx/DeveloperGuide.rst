.. highlight:: rst


Developer Guidelines
====================

The following developer guidelines apply to all code that is to be committed to the `Pele` suite. 

Before adding files for a commit to be pushed to a branch of `PelePhysics`, first run the following formatting tools depending on the type of code:


C++ Code
--------
All C++ code should be processed by the Clang formatter prior to being added for commit.

Run ``clang-format``::

    clang-format -i FILE.cpp
    clang-format -i FILE.H

Will apply all of the correct formatting and make the changes directly to the cpp source files.


Python
------

The tools necessary to format Python code (currently only used within CEPTR) are maintained through Poetry.

1) Run ``ruff format``::

    poetry run ruff format .

This will perform proper formatting of all Python files (including the sorting of imports) to be consistent with all current files.

2) Run ``ruff check``::

    poetry run ruff check .

This will run diagnostics on all the Python files and will list a series of issues that need to be addressed to adhere to current Python best practices. Many of them can be fixed automatically with ``poetry run ruff check --fix .``.

Once all ``ruff`` messages have been addressed, the code will match the `Pele` suite standard.

