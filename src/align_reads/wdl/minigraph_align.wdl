version 1.0

import "gnomix.wdl" as gnomix_t

workflow minigraph_align {

    meta {
        description: "TODO"
    }
    parameter_meta {
        # TODO
    }

    input {
        String out_prefix
        Array[File] batch_vcf_files
        File gnomix_model
        String chrom
        File chainfile
        File refpanel
        File refpanel_index
    }

    ### Output files ####
    output {
        File msp_outfile = merge_gnomix.msp_outfile
        File fb_outfile = merge_gnomix.fb_outfile
    }
}

### Example Task
task minigraph {
    input {
        File pangenome_gfa
        File assembly_fa
        String out_prefix
    }

    String timestamp = ""
    String out_gaf = out_prefix + ".gaf"
    String out_log = out_prefix + "_" + timestamp + ".log"

    command <<<
    minigraph 
        -x asm \
        ~{pangenome_gfa} \
        ~{assembly_fa} \
        -t 4 \
            > ~{out_gaf} \
            2> ~{out_log}
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/beagle:latest"
        memory: "20GB"
    }

    output {
       File out_gaf = out_gaf
       File out_log = out_log
    }
}