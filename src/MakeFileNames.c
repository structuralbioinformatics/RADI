#include "raDI.h"

void MakeFileNames(root,ra,cmap_filename,cmap_off_filename,MI_top_filename,DI_top_filename,DI_filename,gnu_filename,gnu_gif)
int   ra;
char *root,*DI_top_filename,*DI_filename,*MI_top_filename,*cmap_filename,*cmap_off_filename,*gnu_filename,*gnu_gif;
{
//Define plot/output names
        sprintf(cmap_filename,"%s_ra%1d_cmap.dat",root,ra);
        sprintf(cmap_off_filename,"%s_ra%1d_cmapOff.dat",root,ra);
        sprintf(MI_top_filename,"%s_ra%1d_MI_top.dat",root,ra);
        sprintf(DI_top_filename,"%s_ra%1d_DI_top.dat",root,ra);
        sprintf(DI_filename,"%s_ra%1d_DI.out",root,ra );
        sprintf(gnu_filename,"%s_ra%1d_gnu.dat",root,ra );
        sprintf(gnu_gif,"%s_ra%1d_cmap.gif",root,ra );

}
