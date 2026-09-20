# Push misoR to a private GitHub repository

Run these commands inside the package directory, not the parent research directory. They publish only this standalone package. No remote repository has been created by the package setup.

## GitHub website and Git

1. On GitHub, create a repository named `misoR`, choose **Private**, and leave the README, license, and `.gitignore` initialization options unchecked.
2. In a terminal, initialize and commit the package. Replace the directory below if you move it. The Git identity settings are local to this repository.

```sh
cd /Users/trongdatdo/Dat/miso-main/misoR
git init -b main
git config user.name "Dat Do"
git config user.email "dodat.stats@gmail.com"
git add .
git commit -m "Initial misoR package"
```

3. Replace `YOUR_USERNAME` with the repository owner and push:

```sh
git remote add origin https://github.com/YOUR_USERNAME/misoR.git
git push -u origin main
```

Use your configured GitHub credential manager/token for HTTPS authentication, or use an SSH remote if you already have GitHub SSH access.

## Alternative: GitHub CLI

If `gh` is installed, initialize and commit locally as above, then create the private repository and push in one command (skip the website creation and manual remote steps):

```sh
gh auth login
gh repo create misoR --private --source=. --remote=origin --push
```

The [official GitHub CLI documentation](https://cli.github.com/manual/gh_repo_create) describes these options. `--private` selects private visibility, `--source=.` uses this package, and `--push` uploads the commits.

## Subsequent changes

```sh
git add .
git commit -m "Describe your changes"
git push
```

The repository ignores R session files, `.Renviron`, build archives, and check directories. The package author and maintainer are recorded in `DESCRIPTION`. Choose an appropriate distribution license before a future public release; the current notice retains rights for private research development.
