#!/usr/bin/env nextflow
/**
RNA outliers to interactive report workflow
**/

process ResultsToHtml {
    time '1h'
    memory '4 GB'
    cpus 1

    publishDir "$params.output/report/", mode: 'copy'

    input:
        tuple path(samplesheet), path(files)
    output:
        path "*.html"

    script:
        """
        ${CMD_REPORT} bash -c '

        # Determine if samplesheet contains a "genePanel" column
        header=\$(head -n1 "$samplesheet")
        has_gene_panel=false
        gene_panel_index=-1

        IFS=\$"\t" read -r -a header_cols <<< "\$header"
        for i in "\${!header_cols[@]}"; do
            if [[ "\${header_cols[\$i]}" == "genePanel" ]]; then
                has_gene_panel=true
                gene_panel_index=\$i
                echo "Detected genePanel column at index \$gene_panel_index"
            fi
        done

        # Process rows
        tail -n +2 "$samplesheet" | while IFS=\$"\t" read -r -a row; do
            sampleid="\${row[0]}"

            echo "Now processing \${sampleid}:"

            outrider_path="\${sampleid}_result_table_outrider.tsv"
            fraser_path="\${sampleid}_result_table_fraser.tsv"
            mae_path="\${sampleid}_result_mae.tsv"
            output_path="\${sampleid}_report.html"

            # Samples without a bam (OUTRIDER only) have no FRASER results: no report
            if [[ ! -f "\$fraser_path" ]]; then
                echo "No FRASER results for \${sampleid}, skipping report"
                continue
            fi

            # Handle optional MAE file
            if [[ -f "\$mae_path" ]]; then
                echo "\$mae_path"
            else
                mae_path="-"
            fi
            gene_panel_path=""
            
            # Detect genePanel file if column exists
            if [[ "\$has_gene_panel" = true ]]; then
                gene_panel_value="\${row[\$gene_panel_index]}"

                gene_panel_path="\$(basename \$gene_panel_value)"

                # Only check if value is non-empty
                if [[ -n "\$gene_panel_value" ]]; then

                    if [[ -f "\$gene_panel_path" ]]; then
                        echo "Found genePanel file: \$gene_panel_path"
                    else
                        echo "No .txt or .bed genePanel file found for value: \$gene_panel_path"
                        gene_panel_path=""

                    fi
                fi
            fi

            # Build command
            cmd="python3 ${params.report.embedScript} -or \${outrider_path} -fr \${fraser_path} -ma \${mae_path} -t ${params.report.htmlTemplate} -s \${sampleid} -o \${output_path}"

            # Add genePanel argument only when file exists
            if [[ -n "\$gene_panel_path" ]]; then
                cmd="\$cmd -g \$gene_panel_path"
            fi

            echo "Running: \$cmd"
            eval \$cmd

        done
        '
        """
        
}