
import java.nio.file.*


domain   = params.domain
ncbiPath = 'ftp://ftp.ncbi.nlm.nih.gov/genomes'


Channel.fromPath( "${ncbiPath}/refseq/${domain}/assembly_summary.txt" )
    .splitCsv( skip: 2, sep: '\t'  )
    .filter( { (it[11].toLowerCase() == 'complete genome')  ||  (it[4].toLowerCase() == 'representative genome')  ||  (it[4].toLowerCase() == 'reference genome') } )
    .map( {
        def species = it[7]
        def strain  = it[8] - 'strain='
	def status  = it[11].split(' ')[0].toLowerCase()
        if( species.contains( strain ) )
            return [ it[0], it[5], species, status, it[19] - 'ftp://ftp.ncbi.nlm.nih.gov/genomes/' ]
        else
            return [ it[0], it[5], "${species} ${strain}", status, it[19] - 'ftp://ftp.ncbi.nlm.nih.gov/genomes/' ]
    } )
    .set { validGenomes }


process sketch {

    maxForks 4
    maxRetries 3
    tag { "${acc} - ${orgName}" }

    input:
    set val(acc), val(taxId), val(orgName), val(status), val(path) from validGenomes

    output:
    set val(acc), val(taxId), val(status), val(orgName) into dbEntries
    file("${acc}.msh") into outMash
    file("${acc}.fna") into outFasta

    publishDir path: "./${domain}/", mode: 'link'

    script:
    """
    wget -O ${acc}.gz ${ncbiPath}/${path}/${path.split('/').last()}_genomic.fna.gz
    gunzip ${acc}.gz
    $REFERENCE_SEEKER_HOME/share/mash/mash sketch -k 32 -s 10000 ${acc}
    mv ${acc} ${acc}.fna
    """
}


dbEntries.map { "${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}" }
    .collectFile( name: 'db.tsv', storeDir: "./${domain}/", newLine: true )
