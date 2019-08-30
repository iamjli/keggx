# KEGG download

Find ID of organism
`curl -s  "http://rest.kegg.jp/list/organism" | grep "Arabidopsis thaliana"`

Download pathway xml files
`curl "http://rest.kegg.jp/list/pathway/T01001" | cut -f 1 | while read A; do  curl -o "human_KGML/${A:5}.xml" "http://rest.kegg.jp/get/${A}/kgml" ; done`

Download compound ID mapping
`curl "http://rest.kegg.jp/list/compound" > KEGG_compound_ids.txt`


# Reactome download

Download functional interactions (FIs)
`wget http://reactomews.oicr.on.ca:8080/caBigR3WebApp2016/FIsInGene_022717_with_annotations.txt.zip`