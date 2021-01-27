# Mars500 Oral Microbiome [![made-with-Markdown](https://img.shields.io/badge/Made%20with-Markdown-1f425f.svg)](http://commonmark.org)

<img alt="Logo of the project «Mars-500».svg" src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/e3/Logo_of_the_project_%C2%ABMars-500%C2%BB.svg/170px-Logo_of_the_project_%C2%ABMars-500%C2%BB.svg.png" width="170" height="119" align="right">

This repository contains all the codes used to generate files and table reported in the work:

Bacci et al. (2020) [Defining the resilience of the human salivary microbiota by a 520 days longitudinal study in confined environment: the Mars500 mission](https://www.biorxiv.org/content/10.1101/2020.04.08.031401v2) *bioRxiv*

If you use code from this repo, please cite our paper as follows:

```
@article {Bacci2020.04.08.031401,
	author = {Bacci, Giovanni and Mengoni, Alessio and Emiliani, Giovanni and Chiellini, Carolina and Cipriani, Edoardo Giovanni and Bianconi, Giovanna and Canganella, Francesco and Fani, Renato},
	title = {Defining the resilience of the human salivary microbiota by a 520 days longitudinal study in confined environment: the Mars500 mission},
	elocation-id = {2020.04.08.031401},
	year = {2020},
	doi = {10.1101/2020.04.08.031401},
	publisher = {Cold Spring Harbor Laboratory},
	eprint = {https://www.biorxiv.org/content/early/2020/05/12/2020.04.08.031401.full.pdf},
	journal = {bioRxiv}
}
```

This work is part of the Mars500 mission, the first prolonged isolation experiment involving human volunteers conducted between 2007 and 2011 by three space agencies of different countries: Russia (Roscosmos), Europe (ESA), and China (CNSA). The mission was part of the European Programme for Life and Physical Sciences (ELIPS) with the aim of preparing astronauts for future missions to the Moon and Mars.

Details about the project can be found [here](https://www.esa.int/Science_Exploration/Human_and_Robotic_Exploration/Mars500/Mars500_study_overview)

## Project abstract

The human microbiota plays several roles in health and disease but is often difficult to determine which part is in intimate relationships with the host _vs._ the occasional presence. During the Mars500 mission, six crewmembers lived completely isolated from the outer world for 520 days following standardized diet regimes. The mission constitutes the first spaceflight simulation to Mars and was a unique experiment to determine, in a longitudinal study design, the composition and importance of the resident _vs._ a more variable microbiota--the fraction of the human microbiota that changes in time and according to environmental conditions--in humans. Here we report the characterization of the salivary microbiota from 88 samples taken during and after the mission for a total of 720 days. Amplicon sequencing of the V3-V4 region of 16S rRNA gene was performed and results were analyzed monitoring the diversity of the microbiota while evaluating the effect of the three main variables present in the experimental system: time, diet, and individuality of each subject. Results, though showing statistically significant effects of all three variables, highlighted a main contribution of salivary microbiota personalized features, that is an individual-based resilience of the microbiota. Such findings open the way to consider salivary microbiota under the light of a pronounced personalization even after sharing the same physical space for more than a year.

## File description

The two main folders of the repository are called `SM` and `main` and refer to supplementary materials and main text, respectively. The content of the folders is organised as follows:

* `img` folder: images generated during `knitr` rendering are saved into this folder in `png` format

* `data` folder: contains data that is loaded into the R environment to generate figures and external scripts used for time-consuming tasks

* `tab` folder (SM only): supplementary tables are saved into this folder during `knitr` rendering of supplementary materials. Tables are embedded in the main text and are not saved when rendering the document.

* `main/manuscript.Rmd`: R Markdown file for generating figures and table reported in the main manuscript

* `SM/supplementary.Rmd`: R Markdown file for generating figures and table reported in supplementary materials

## Instructions

R Markdown files can be rendered into `html` by using the `knitr` plugin integrated into [RStudio](https://rstudio.com/?_ga=2.50552553.1339302526.1611745574-1183453795.1578408315). Alternatively, one can directly render documents from terminal by launching this command:

```shell
R -e "rmarkdown::render('<RMarkdown_file>',output_file='<output_file>')"
```

by replacing `<RMarkdown_file>` with the name of the file to render and `<output_file>` with the name of the output file that will be generated.

It is important to note that the `html` documents are generated using `pandoc` package integrated with RStudio. In Debian/Linux systems the library can be found in `/usr/lib/rstudio/bin/pandoc/` and could differ from the version installed at system level. If, during the rendering of the document with `knitr`, you get a message like:

```shell
Error: pandoc version 1.12.3 or higher is required and was not found (see the help page ?rmarkdown::pandoc_available).
```

then you probably need to tell `R` where the `pandoc` library integrated with RStudio is located in your system. You can do this with the command `Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/bin/pandoc/')` or by using the two build scripts provided together with the markdown files.
