library(data.table)

wdir = "../sim/"
odir = "../doc/"
fdir = "../fig/"
setwd(wdir)

qfolderlist = c('3_0.05_0.05/', '2_0.05_0.05/', '1_0.05_0.05/')
Locilist = c("Locus1/", "Locus2/", "Locus3/")

ngem_pval_lst = list()
ncojo_pval_lst = list()
nclump_pval_lst = list()
nsharepro_pval_lst = list()
for (folder in qfolderlist) {
    for (rt in c('CL','CM','CS')) {
        gxegem = c()
        gxecojo = c()
        gxeclump = c()
        gxesharepro = c()
        for (lo in Locilist) {
            for (x in 1:50) {
                gem = read.table(paste0(wdir,lo,folder,rt,x), header=T)
                rownames(gem) = gem$SNPID
                gem$cojop = 0
                gem$clumpp = 0
                gem$shp = 0
                if (paste0(rt,x,'.jma.cojo') %in% list.files(paste0(wdir,lo,folder))) {
                    cojo = fread(paste0(wdir,lo,folder,rt,x,'.jma.cojo'))
                    gem[cojo$SNP,'cojop'] = 1
                }
                if (paste0(rt,x,'.clumped') %in% list.files(paste0(wdir,lo,folder))) {
                    clump = fread(paste0(wdir,lo,folder,rt,x,'.clumped'))
                    gem[clump$SNP,'clumpp'] = 1
                }
                SH = fread(paste0(wdir,lo,folder,rt,x,'.sharepro.gxe.txt'))
                gem[SH$cs!=0,'shp'] = 1
                gxegem = c(gxegem, gem$P_Value_Interaction)
                gxecojo = c(gxecojo, gem$cojop)
                gxeclump = c(gxeclump, gem$clumpp)
                gxesharepro = c(gxesharepro, gem$shp)
            }
        }
        name = paste0(folder, rt)
        ngem_pval_lst[[name]] = gxegem
        ncojo_pval_lst[[name]] = gxecojo
        nclump_pval_lst[[name]] = gxeclump
        nsharepro_pval_lst[[name]] = gxesharepro
    }
}

save(ngem_pval_lst, ncojo_pval_lst, nclump_pval_lst, nsharepro_pval_lst, file='../doc/calibration.Rdata')