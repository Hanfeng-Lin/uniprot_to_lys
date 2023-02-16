from lxml import etree
import requests
import pandas as pd


def uniprot_to_pdb(uniprot_id):
    # Get the XML content from the URL
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.xml"
    response = requests.get(url)
    xml_string = response.content

    # Parse the XML content
    parsed_xml = etree.fromstring(xml_string)  # parsed_xml = etree.parse("P29466.xml")
    namespaces = {'ns': 'http://uniprot.org/uniprot'}

    # Find the chains value for the given PDB ID
    # X-ray only
    seq = parsed_xml.xpath("//ns:sequence/text()", namespaces=namespaces)[0]
    name = parsed_xml.xpath("//ns:gene/ns:name[@type='primary']/text()", namespaces=namespaces)[0]
    domain_dict = {}
    domain_features = parsed_xml.xpath("//ns:feature[@type='domain' or @type='short sequence motif']",
                                       namespaces=namespaces)
    df1 = pd.DataFrame()

    for feature in domain_features:
        description = feature.get("description")
        domain_begin = int(feature.xpath(".//ns:location/ns:begin/@position", namespaces=namespaces)[0])
        domain_end = int(feature.xpath(".//ns:location/ns:end/@position", namespaces=namespaces)[0])
        domain_dict[description] = [domain_begin, domain_end]

    print(domain_dict)

    pdb_features = parsed_xml.xpath("//ns:entry/ns:dbReference[@type='PDB']", namespaces=namespaces)
    if pdb_features:
        for pdb in pdb_features:
            method = pdb.xpath("./ns:property[@type='method']/@value", namespaces=namespaces)[0]
            pdb_id = pdb.get("id")
            if method == "X-ray":
                chains = pdb.xpath("./ns:property[@type='chains']/@value", namespaces=namespaces)[0]
                for chain in chains.split(","):
                    [pdb_begin, pdb_end] = list(map(int, chain.split("=")[1].split("-")))
                    pdb_seq = seq[pdb_begin - 1:pdb_end]
                    lys_number = pdb_seq.count("K")
                    # Find what domain does this structure contain
                    if domain_dict:
                        for domain, [domain_start, domain_end] in domain_dict.items():
                            if all(elem in range(pdb_begin, pdb_end + 1) for elem in
                                   range(domain_start, domain_end + 1)):  # check domain_range is a subset of pdb_range
                                print(
                                    f"{uniprot_id} {name} Chains for PDB ID {pdb_id}:  {chain} chain contains lys: {lys_number}, contains domain:{domain} {[domain_start, domain_end]}  chain sequence: {pdb_seq}")
                                df2 = pd.DataFrame(
                                    [uniprot_id, name, pdb_id, chain, domain, domain_start, domain_end, lys_number,
                                     pdb_seq]).transpose()
                                df1 = pd.concat([df1, df2])
                    else:
                        df2 = pd.DataFrame(
                            [uniprot_id, name, pdb_id, chain, "N/A", "", "", lys_number,
                             pdb_seq]).transpose()
                        df1 = pd.concat([df1, df2])
            else:
                print(pdb_id, method)
    else:
        df1 = pd.DataFrame([uniprot_id, name, "N/A", "", "", "", "", "", ""]).transpose()

    return df1


uniprot_list = ["Q2M2I8", "P00519", "P42684", "P35348", "P31749", "P31751", "Q9Y243", "Q9UM73", "P37840", "P10275",
                "O14965", "Q96GD4", "P30530", "O14874", "P10415", "P41182", "Q07817", "P07437", "P51451", "Q9NSY1",
                "P36894", "P15056", "P25440", "Q15059", "O60885", "Q9NPI1", "Q9H8M2", "Q06187", "O43683", "O60566",
                "Q8N5S9", "Q92793", "P51686", "P35613", "Q12834", "O00311", "P06493", "Q15131", "Q9UQ88", "P21127",
                "Q9NYV4", "Q14004", "O94921", "Q00536", "Q00537", "Q07002", "Q9BWU1", "P24941", "P11802", "Q00535",
                "Q00534", "P50613", "P49336", "P50750", "O14757", "O14578", "P10721", "P49759", "P08581", "Q8NI60",
                "Q96D53", "P29762", "P29373", "Q96SW2", "P41240", "P48729", "P48730", "P49674", "Q16678", "P53355",
                "Q16832", "O75530", "O00418", "P00533", "P19525", "Q9P2K8", "P06730", "Q09472", "P21709", "P29317",
                "P29320", "P29323", "P54753", "P54760", "O15197", "P03372", "O75460", "P11474", "Q15910", "Q05397",
                "P16591", "P11362", "P21802", "P62942", "P17948", "P36888", "P06241", "O14976", "Q92830", "P49840",
                "P49841", "Q8TF76", "Q13547", "Q92769", "O15379", "Q9UBN7", "Q9BY41", "P04626", "Q86Z02", "P04035",
                "O60760", "Q00613", "P07900", "P42858", "P14902", "P08069", "P51617", "Q9Y616", "Q9NWZ3", "Q08881",
                "P23458", "O60674", "P52333", "P01116", "O95835", "P06239", "P53667", "P53671", "Q5S007", "Q13133",
                "P55055", "P07948", "Q13163", "Q13233", "Q16584", "Q12852", "Q9NYL2", "Q5TCX8", "O43318", "Q92918",
                "Q12851", "Q8IVH8", "Q9Y4K4", "P53779", "P53778", "Q16659", "Q13164", "P45983", "P45984", "P49137",
                "Q16644", "Q8IW41", "Q7KZI7", "P27448", "Q96L34", "O60307", "Q07820", "Q00987", "Q02750", "P36507",
                "Q14680", "Q12866", "P14174", "Q8N4C8", "Q9HBH9", "Q03111", "P01106", "Q96PY6", "P51955", "P51956",
                "P51957", "Q8TD19", "Q9UBE8", "O60285", "Q16539", "Q15759", "O15264", "Q13153", "O96013", "P09874",
                "Q460N5", "Q9UGN5", "Q9Y6F1", "Q86U86", "Q92831", "O43924", "Q8N165", "Q15118", "Q15119", "Q15120",
                "Q9NZQ7", "P42336", "P48736", "Q9P1W9", "P48426", "P78356", "Q8TBX8", "Q99640", "Q6P5Z2", "P53350",
                "O00444", "Q86YV5", "Q13131", "P54646", "P41743", "O14744", "Q14289", "Q13882", "Q9BVS4", "Q13546",
                "O43353", "Q16186", "Q15418", "P51812", "O75676", "Q9UK32", "P23443", "Q96S38", "Q52WX2", "O75533",
                "Q96BR1", "Q06124", "Q9H0K1", "Q9Y2K2", "Q8IXJ6", "P19634", "Q9UBY0", "Q6AI14", "Q96T83", "Q4ZJI4",
                "P51531", "P51532", "Q9NRH2", "Q07889", "P12931", "P40763", "O94804", "Q9UEE5", "O94768", "Q86UX6",
                "Q9BYT3", "Q8TDR2", "Q15208", "Q13043", "Q8N2I9", "Q15022", "Q9Y6A5", "Q9UL54", "Q9H2K8", "P10636",
                "Q8TEA7", "Q9UHD2", "P42680", "Q96S53", "P36897", "Q13470", "Q07912", "Q96RU7", "O15164", "P04629",
                "Q16620", "Q16288", "Q96QT4", "P33981", "P29597", "P14679", "Q8TAS1", "O75385", "Q6PHR2", "P35968",
                "P40337", "P61964", "P30291", "P07947"]

df1 = pd.DataFrame()
for uniprot in uniprot_list:
    df2 = uniprot_to_pdb(uniprot)
    df1 = pd.concat([df1, df2])

df1.to_excel("uniprot_to_lys.xlsx", index=False)
