<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- [Developer Guide](#developer-guide)
  - [Additional requirements](#additional-requirements)
  - [Downloading MODA](#downloading-moda)
  - [Installing Git hooks](#installing-git-hooks)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Developer Guide 

This guide is aimed at developers wishing to modify or contribute to the program.

## Additional requirements

In addition to the [requirements](/README.md#requirements) listed in the User Guide, you'll also need to [install Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

## Downloading MODA

- Open a terminal in a desired folder and run `git clone https://github.com/luphysics/MODA.git`.
- The code will download as a folder named `MODA`.

## Installing Git hooks

Git hooks are used to automatically perform tasks when a commit is made. MODA uses `doctoc` to add the table of contents to markdown files.

Commit your current work, if there are changes. Then open a terminal in the `MODA` folder and run:

```
pip install pre-commit --user   # Installs the pre-commit tool.
python -m pre-commit install    # Adds the Git hooks to the repository.
```

On Windows, also run `git config core.safecrlf false` in the `MODA` folder. This prevents a circular problem where Git cannot commit because it converts line endings to CRLF but `doctoc` converts line endings back to LF.

Now that the Git hooks are installed, they will automatically run every time a commit changes relevant files.

> :warning: When a pre-commit hook changes files, you'll need to add files and commit again.