function myprint()
{
	for i in `cat /afs/cern.ch/work/f/fnemes/main_workspace_github_ssh_4/Projects/TOTEM_Projects/Physics_projects/Physics_analysis/Proton_proton_scattering/Elastic_scattering/Projects/2023/E_CM_900_GeV_beta_star_100_m/General_settings/List_of_runs.txt` ; do echo "\begin{frame}\tiny" ; echo $i ; grep $1.*$i status.tex ; echo "\end{frame}" ; echo  ; done
}


myprint Left_bottom.
myprint Left_top
