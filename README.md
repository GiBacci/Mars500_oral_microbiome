# Mars500 Oral Microbiome [![made-with-Markdown](https://img.shields.io/badge/Made%20with-Markdown-1f425f.svg)](http://commonmark.org)

<img alt="Logo of the project «Mars-500».svg" src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/e3/Logo_of_the_project_%C2%ABMars-500%C2%BB.svg/170px-Logo_of_the_project_%C2%ABMars-500%C2%BB.svg.png" width="170" height="119" align="right">

This repository contains all the codes used to generate files and table reported in the work:

Bacci et al. (2020) [Defining the resilience of the human salivary microbiota by a 520 days longitudinal study in confined environment: the Mars500 mission](https://www.biorxiv.org/content/10.1101/2020.04.08.031401v2) *bioRxiv*

The project is part of the Mars500 mission, the first prolonged isolation experiment involving human volunteers conducted between 2007 and 2011 by three space agencies of different countries: Russia (Roscosmos), Europe (ESA), and China (CNSA). The mission was part of the European Programme for Life and Physical Sciences (ELIPS) with the aim of preparing astronauts for future missions to the Moon and Mars.

Details about the project can be found [here](https://www.esa.int/Science_Exploration/Human_and_Robotic_Exploration/Mars500/Mars500_study_overview)

## Project abstract

The human microbiota plays several roles in health and disease but is often difficult to determine which part is in intimate relationships with the host _vs._ the occasional presence. During the Mars500 mission, six crewmembers lived completely isolated from the outer world for 520 days following standardized diet regimes. The mission constitutes the first spaceflight simulation to Mars and was a unique experiment to determine, in a longitudinal study design, the composition and importance of the resident _vs._ a more variable microbiota---the fraction of the human microbiota that changes in time and according to environmental conditions---in humans. Here we report the characterization of the salivary microbiota from 88 samples taken during and after the mission for a total of 720 days. Amplicon sequencing of the V3-V4 region of 16S rRNA gene was performed and results were analyzed monitoring the diversity of the microbiota while evaluating the effect of the three main variables present in the experimental system: time, diet, and individuality of each subject. Results, though showing statistically significant effects of all three variables, highlighted a main contribution of salivary microbiota personalized features, that is an individual-based resilience of the microbiota. Such findings open the way to consider salivary microbiota under the light of a pronounced personalization even after sharing the same physical space for more than a year.

## File description

The two main fodlers of the repository are called `SM` and `main` and refer to supplementary materials and main text, respectively. The content of the foldes is organised as follows:

1. `img` folder: images generated during `knitr` rendering are saved into this folder in `png` format.

2. `data` folder: contains data that is loaded into the R environment to generate figures and external scripts used for time-consuming tasks
    
    <details>
     <summary><b>Details</b></summary>
    1. `SM/data/cp` folder: contains the code used for change point analysis<br>
    2. `SM/data/taxa` foder: contains the code for database comparison<br>
    3. `SM/data/all_seqtab_nochim.rds`: ASV table without filtering<br>
    4. `SM/data/phylo_obj.rds`: phyloseq object containing the final counts, metadat, and taxonoic assignments
    </details>

3. `tab` folder (SM only): supplementary tables are saved into this folder during `knitr` rendering of supplementary materials. Tables are embedded in the main text and are not saved when rendering the document.
