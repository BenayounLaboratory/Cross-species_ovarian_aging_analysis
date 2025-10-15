#!/bin/bash

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id YF_1 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/AC \
                 --sample MHK333 \
                 --localmem 64 \
                 --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id YF_2 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/AC \
                 --sample MHK334 \
                 --localmem 64 \
                 --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id OF_1 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/AC \
                 --sample MHK335 \
                 --localmem 64 \
                 --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id OF_2 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/AC \
                 --sample MHK336 \
                 --localmem 64 \
                 --localcores 12
                 
                 
/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id Foxl2_wt \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/Foxl2/10xFQ \
                 --sample MHK318 \
                 --localmem 64 \
                 --localcores 12
                 
/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id Foxl2_wt_young_1 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/Foxl2/10xFQ \
                 --sample MHK337 \
                 --localmem 64 \
                 --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id Foxl2_wt_old_1 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/Foxl2/10xFQ \
                 --sample MHK338 \
                 --localmem 64 \
                 --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id Foxl2_wt_young_2 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/Foxl2/10xFQ \
                 --sample MHK360 \
                 --localmem 64 \
                 --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id Foxl2_wt_old_2 \
                 --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
                 --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/Foxl2/10xFQ \
                 --sample MHK361 \
                 --localmem 64 \
                 --localcores 12
               
/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id CTL_3m_30d_1 \
         --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
         --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/VCD \
         --sample MHK305 \
         --localmem 64 \
         --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id CTL_3m_90d_1 \
         --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
         --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/VCD \
         --sample MHK307 \
         --localmem 64 \
         --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id CTL_10m_30d_1 \
         --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
         --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/VCD \
         --sample MHK309 \
         --localmem 64 \
         --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id CTL_10m_90d_1 \
         --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
         --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/VCD \
         --sample MHK311 \
         --localmem 64 \
         --localcores 12
          
                 
/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id CTL_3m_30d_2 \
         --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
         --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/VCD \
         --sample MHK320 \
         --localmem 64 \
         --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id CTL_3m_90d_2 \
         --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
         --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/VCD \
         --sample MHK328 \
         --localmem 64 \
         --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id CTL_10m_30d_2 \
         --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
         --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/VCD \
         --sample MHK322 \
         --localmem 64 \
         --localcores 12

/home/minhooki/Softwares/cellranger-7.1.0/cellranger count --id CTL_10m_90d_2 \
         --transcriptome /home/benayoun/Softwares/cellranger-6.0.2/Reference/refdata-gex-mm10-2020-A \
         --fastqs /mnt/CZI_data_1/Data/Benayoun_lab/1_Rawdata/VCD \
         --sample MHK330 \
         --localmem 64 \
         --localcores 12