# How to publish a new release

When a sufficient number of changes have been accumulated, a new release
can be made. There are three places in which a release must be published:
github, pypi and bioconda.

**Note:** `panfeed` follow's the [semantic versioning](https://semver.org/) convention.

## 1) Add the release version to the repository

In the `main` branch, 
edit the `panfeed/__init__.py` file to indicate the target release version (e.g. `X.X.X`).
Then do `git add panfeed/__init__.py` followed by `git commit -m "Version bump"`.

## 2) Build panfeed for pypi and github

This step requires the `hatchling` python package, and that the `panfeed` directory
does not contain any file that is not part of the git repository; doing this on a fresh
repository is the easiest way to avoid adding unwanted files.

    python3 -m build

The above command should create two files in the `dist` directory:
`panfeed-X.X.X-py3-none-any.whl` and `panfeed-X.X.X.tar.gz`.

## 3) Make a github release

Apply a tag to identify the release in the git history by doing: `git tag X.X.X`
(where `X.X.X` is the target version), followed by `git push --tags`.
Go to [panfeed's release page](https://github.com/microbial-pangenomes-lab/panfeed/releases)
, then click on "tags", click on the tag you just pushed, and finally click on
"Make new release from tag".
Fill the various fields following what has been done for previous releases
(e.g. [this one](https://github.com/microbial-pangenomes-lab/panfeed/releases/tag/1.5.0)),
and add the two files generated by the previous step, then publish the new release.

## 4) Upload on pypi

This step requires the `twine` python package, an account on [pypi.org](https://pypi.org/) and
that the account has been added as a maintainer for [panfeed](https://pypi.org/project/panfeed/).

    python3 -m twine upload dist/*

## 5) Upload on bioconda

**NOTE:** it might be that bioconda will automatically create a pull request when `panfeed` is updated
on pypi. So please wait up to a day before opening a pull request manually.

This step requires that the upload on pypi was successful, as it relies
on the new release to be present there. It also assumes you have a fork of the huge
[`bioconda-recipes`](https://github.com/bioconda/bioconda-recipes/) repository and at least ~8GB
of RAM and a linux system. [Follow the steps here](https://bioconda.github.io/contributor/setup.html)
to set up your own fork.

Set up a bioconda build environment by creating a `bioconda` environment in your local `conda` installation:

    mamba create -n bioconda -c conda-forge -c bioconda bioconda-utils grayskull
    conda activate bioconda

Then `cd` into your local copy of `bioconda-recipes` and update it to the latest upstream commits:

    git checkout master
    git pull upstream master
    git push origin master

Then create a new branch to update the recipe:

    git checkout -b update_panfeed_XXX

Then use `grayskull` to update the recipe:

    cd recipes
    grayskull pypi panfeed

Then use `git add -p` to accept/reject changes. Please just add changes to the version and the
sha256 hash of the package file. The other changes are just cosmetic, unless a new subcommand
or dependency has been added. Create a commit with `git commit`, followed by `git stash` to get
rid of the unwanted changes.

Test the updated recipes with:

    cd ..
    bioconda-utils lint --packages panfeed
    bioconda-utils build --packages panfeed

If everything looks good, do `git push --set-upstream origin update_panfeed_XXX`.
The remote server will likely provide a link to directly open a pull request based on these
changes, otherwise visit your fork of `bioconda-recipes` and there should be a button somewhere to
create a pull request. Follow the indications in the template pull request to provide
information, then wait for the remote tests to complete and add a command saying:
`@bioconda-bot please add label`. A bioconda maintainer will eventually merge and the new version will
then appear on bioconda.

## 6) Bump panfeed's repository

On `panfeed`'s `main` branch, edit the `panfeed/__init__.py` so that it's clear the next version
is a draft (e.g. if latest release is `1.5.0` you could do `1.5.1-dev`), followed by
`git add panfeed/__init__.py`, `git commit -m "Dev bump"` and `git push`.
