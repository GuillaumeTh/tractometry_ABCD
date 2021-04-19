#!/usr/bin/env nextflow

log.info "Tractometer 2.0 pipeline"
log.info "========================"
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "Input: $params.input"
root = file(params.input)
/* Watch out, files are ordered alphabetically in channel */

metrics = Channel
        .fromFilePairs("$root/**/*{ad,fa,fodf,md,nufo,rd}.nii.gz",
                       size: 6,
                       maxDepth:1,
                       flat: true) {it.parent.name}

(fodf_for_afd, metrics_for_along) = metrics
    .map{sid, ad, fa, fodf, md, nufo, rd -> [tuple(sid, fodf),
                                        tuple(sid, ad, fa, md, rd, nufo)]}
    .separate(2)

bundles = Channel
        .fromFilePairs("$root/**/*.trk",
                       size: -1,
                       maxDepth:1,
                       flat: true) {it.parent.name}

bundles
    .map{[it[0], it[1..-1]]}
    .into{bundles_for_afd;bundles_for_along}

fodf_for_afd
    .join(bundles_for_afd)
    .set{fodf_bundles_for_afd}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

process Compute_AFD{

input:
    set sid, file(fodf), file(bundles) from fodf_bundles_for_afd

output:
    set sid, "*afd.nii.gz" into afd_out

shell:
'''
for i in !{bundles};
do
scil_compute_mean_fixel_afd_from_bundles.py ${i} !{fodf} !{sid}_${i}__afd.nii.gz !{sid}_${i}__radfodf.nii.gz
done
'''
}

metrics_for_along
    .join(afd_out)
    .join(bundles_for_along)
    .set{metrics_afd_for_along}

process Compute_Mean_Bundles{
echo true
input:
    set sid, file(ad), file(fa), file(md), file(rd), file(nufo), file(afd), file(bundles) from metrics_afd_for_along

output:
    set val("out"), "${sid}__all_metrics_along.json" into metrics_per_subject

shell:
'''
ad=!{ad}
mv !{ad} ${ad/!{sid}__/}
fa=!{fa}
mv !{fa} ${fa/!{sid}__/}
md=!{md}
mv !{md} ${md/!{sid}__/}
rd=!{rd}
mv !{rd} ${rd/!{sid}__/}
nufo=!{nufo}
mv !{nufo} ${nufo/!{sid}__/}
for i in !{bundles};
do
mv ${i} ${i/!{sid}__/}
cur_bdl=${i/!{sid}__/}
mv !{sid}_${i}__afd.nii.gz ${cur_bdl/.trk/}_afd.nii.gz
scil_compute_bundle_mean_std.py ${i/!{sid}__/} ${ad/!{sid}__/} ${fa/!{sid}__/} ${md/!{sid}__/} ${rd/!{sid}__/} ${nufo/!{sid}__/} ${cur_bdl/.trk/}_afd.nii.gz > !{sid}_${cur_bdl/.trk/}__mean_metrics_along.json
done
scil_merge_json.py *mean_metrics_along.json !{sid}__all_metrics_along.json --add_parent_key !{sid}
'''
}

metrics_per_subject
    .groupTuple()
    .set{all_jsons}

process Prepare_output{
input:
    set sid, file(jsons) from all_jsons

output:
    file "subjects_merged.csv"

shell:
'''
scil_merge_json.py !{jsons} subjects_merged.json
convert_to_csv.py subjects_merged.json subjects_merged.csv
'''
}
