// Map: ucsc name to links:
def get_genome_ensembl() {
	def list_genome_ensembl = [
"acanthochromis_polyacanthus","accipiter_nisus","ailuropoda_melanoleuca","amazona_collaria","amphilophus_citrinellus","amphiprion_ocellaris","amphiprion_percula","anabas_testudineus","anas_platyrhynchos","anas_platyrhynchos_platyrhynchos","anas_zonorhyncha","ancestral_alleles","anolis_carolinensis","anser_brachyrhynchus","anser_cygnoides","aotus_nancymaae","apteryx_haastii","apteryx_owenii","apteryx_rowi","aquila_chrysaetos_chrysaetos","astatotilapia_calliptera","astyanax_mexicanus","astyanax_mexicanus_pachon","athene_cunicularia","balaenoptera_musculus","betta_splendens","bison_bison_bison","bos_grunniens","bos_indicus_hybrid","bos_mutus","bos_taurus","bos_taurus_hybrid","bubo_bubo","buteo_japonicus","caenorhabditis_elegans","cairina_moschata_domestica","calidris_pugnax","calidris_pygmaea","callithrix_jacchus","callorhinchus_milii","camarhynchus_parvulus","camelus_dromedarius","canis_lupus_dingo","canis_lupus_familiaris","canis_lupus_familiarisbasenji","canis_lupus_familiarisgreatdane","capra_hircus","capra_hircus_blackbengal","carassius_auratus","carlito_syrichta","castor_canadensis","catagonus_wagneri","catharus_ustulatus","cavia_aperea","cavia_porcellus","cebus_capucinus","cercocebus_atys","cervus_hanglu_yarkandensis","chelonoidis_abingdonii","chelydra_serpentina","chinchilla_lanigera","chlorocebus_sabaeus","choloepus_hoffmanni","chrysemys_picta_bellii","chrysolophus_pictus","ciona_intestinalis","ciona_savignyi","clupea_harengus","colobus_angolensis_palliatus","corvus_moneduloides","cottoperca_gobio","coturnix_japonica","cricetulus_griseus_chok1gshd","cricetulus_griseus_crigri","cricetulus_griseus_picr","crocodylus_porosus","cyanistes_caeruleus","cyclopterus_lumpus","cynoglossus_semilaevis","cyprinodon_variegatus","cyprinus_carpio","cyprinus_carpio_germanmirror","cyprinus_carpio_hebaored","cyprinus_carpio_huanghe","danio_rerio","dasypus_novemcinctus","delphinapterus_leucas","denticeps_clupeoides","dicentrarchus_labrax","dipodomys_ordii","dromaius_novaehollandiae","drosophila_melanogaster","echeneis_naucrates","echinops_telfairi","electrophorus_electricus","eptatretus_burgeri","equus_asinus_asinus","equus_caballus","erinaceus_europaeus","erpetoichthys_calabaricus","erythrura_gouldiae","esox_lucius","falco_tinnunculus","felis_catus","ficedula_albicollis","fukomys_damarensis","fundulus_heteroclitus","gadus_morhua","gallus_gallus","gambusia_affinis","gasterosteus_aculeatus","geospiza_fortis","gopherus_agassizii","gopherus_evgoodei","gorilla_gorilla","gouania_willdenowi","haplochromis_burtoni","heterocephalus_glaber_female","heterocephalus_glaber_male","hippocampus_comes","homo_sapiens","hucho_hucho","ictalurus_punctatus","ictidomys_tridecemlineatus","jaculus_jaculus","junco_hyemalis","kryptolebias_marmoratus","labrus_bergylta","larimichthys_crocea","lates_calcarifer","laticauda_laticaudata","latimeria_chalumnae","lepidothrix_coronata","lepisosteus_oculatus","leptobrachium_leishanense","lonchura_striata_domestica","loxodonta_africana","lynx_canadensis","macaca_fascicularis","macaca_mulatta","macaca_nemestrina","malurus_cyaneus_samueli","manacus_vitellinus","mandrillus_leucophaeus","marmota_marmota_marmota","mastacembelus_armatus","maylandia_zebra","meleagris_gallopavo","melopsittacus_undulatus","meriones_unguiculatus","mesocricetus_auratus","microcebus_murinus","microtus_ochrogaster","mola_mola","monodelphis_domestica","monodon_monoceros","monopterus_albus","moschus_moschiferus","mus_caroli","mus_musculus","mus_musculus_129s1svimj","mus_musculus_aj","mus_musculus_akrj","mus_musculus_balbcj","mus_musculus_c3hhej","mus_musculus_c57bl6nj","mus_musculus_casteij","mus_musculus_cbaj","mus_musculus_dba2j","mus_musculus_fvbnj","mus_musculus_lpj","mus_musculus_nodshiltj","mus_musculus_nzohlltj","mus_musculus_pwkphj","mus_musculus_wsbeij","mus_pahari","mus_spicilegus","mus_spretus","mustela_putorius_furo","myotis_lucifugus","myripristis_murdjan","naja_naja","nannospalax_galili","neogobius_melanostomus","neolamprologus_brichardi","neovison_vison","nomascus_leucogenys","notamacropus_eugenii","notechis_scutatus","nothobranchius_furzeri","nothoprocta_perdicaria","numida_meleagris","ochotona_princeps","octodon_degus","oncorhynchus_kisutch","oncorhynchus_mykiss","oncorhynchus_tshawytscha","oreochromis_aureus","oreochromis_niloticus","ornithorhynchus_anatinus","oryctolagus_cuniculus","oryzias_javanicus","oryzias_latipes","oryzias_latipes_hni","oryzias_latipes_hsok","oryzias_melastigma","oryzias_sinensis","otolemur_garnettii","otus_sunia","ovis_aries","ovis_aries_rambouillet","pan_paniscus","pan_troglodytes","panthera_leo","panthera_pardus","panthera_tigris_altaica","papio_anubis","parambassis_ranga","paramormyrops_kingsleyae","parus_major","pavo_cristatus","pelodiscus_sinensis","pelusios_castaneus","periophthalmus_magnuspinnatus","peromyscus_maniculatus_bairdii","petromyzon_marinus","phascolarctos_cinereus","phasianus_colchicus","phocoena_sinus","physeter_catodon","piliocolobus_tephrosceles","podarcis_muralis","poecilia_formosa","poecilia_latipinna","poecilia_mexicana","poecilia_reticulata","pogona_vitticeps","pongo_abelii","procavia_capensis","prolemur_simus","propithecus_coquereli","pseudonaja_textilis","pteropus_vampyrus","pundamilia_nyererei","pygocentrus_nattereri","rattus_norvegicus","rhinolophus_ferrumequinum","rhinopithecus_bieti","rhinopithecus_roxellana","saccharomyces_cerevisiae","saimiri_boliviensis_boliviensis","salarias_fasciatus","salmo_salar","salmo_trutta","salvator_merianae","sander_lucioperca","sarcophilus_harrisii","sciurus_vulgaris","scleropages_formosus","scophthalmus_maximus","serinus_canaria","seriola_dumerili","seriola_lalandi_dorsalis","sinocyclocheilus_anshuiensis","sinocyclocheilus_grahami","sinocyclocheilus_rhinocerous","sorex_araneus","sparus_aurata","spermophilus_dauricus","sphaeramia_orbicularis","sphenodon_punctatus","stachyris_ruficeps","stegastes_partitus","strigops_habroptila","strix_occidentalis_caurina","struthio_camelus_australis","suricata_suricatta","sus_scrofa","sus_scrofa_bamei","sus_scrofa_berkshire","sus_scrofa_hampshire","sus_scrofa_jinhua","sus_scrofa_landrace","sus_scrofa_largewhite","sus_scrofa_meishan","sus_scrofa_pietrain","sus_scrofa_rongchang","sus_scrofa_tibetan","sus_scrofa_usmarc","sus_scrofa_wuzhishan","taeniopygia_guttata","takifugu_rubripes","terrapene_carolina_triunguis","tetraodon_nigroviridis","theropithecus_gelada","tupaia_belangeri","tursiops_truncatus","urocitellus_parryii","ursus_americanus","ursus_maritimus","ursus_thibetanus_thibetanus","varanus_komodoensis","vicugna_pacos","vombatus_ursinus","vulpes_vulpes","xenopus_tropicalis","xiphophorus_couchianus","xiphophorus_maculatus","zalophus_californianus","zonotrichia_albicollis","zosterops_lateralis_melanops","actinidia_chinensis","aegilops_tauschii","amborella_trichopoda","ananas_comosus","arabidopsis_halleri","arabidopsis_lyrata","arabidopsis_thaliana","arabis_alpina","asparagus_officinalis","beta_vulgaris","brachypodium_distachyon","brassica_napus","brassica_oleracea","brassica_rapa","camelina_sativa","cannabis_sativa_female","capsicum_annuum","chara_braunii","chenopodium_quinoa","chlamydomonas_reinhardtii","chondrus_crispus","citrullus_lanatus","citrus_clementina","coffea_canephora","corchorus_capsularis","cucumis_melo","cucumis_sativus","cyanidioschyzon_merolae","cynara_cardunculus","daucus_carota","dioscorea_rotundata","eragrostis_curvula","eragrostis_tef","eucalyptus_grandis","eutrema_salsugineum","galdieria_sulphuraria","glycine_max","gossypium_raimondii","helianthus_annuus","hordeum_vulgare","hordeum_vulgare_goldenpromise","hordeum_vulgare_tritex","ipomoea_triloba","juglans_regia","kalanchoe_fedtschenkoi","leersia_perrieri","lupinus_angustifolius","malus_domestica_golden","manihot_esculenta","marchantia_polymorpha","medicago_truncatula","musa_acuminata","nicotiana_attenuata","nymphaea_colorata","olea_europaea_sylvestris","oryza_barthii","oryza_brachyantha","oryza_glaberrima","oryza_glumipatula","oryza_indica","oryza_longistaminata","oryza_meridionalis","oryza_nivara","oryza_punctata","oryza_rufipogon","oryza_sativa","ostreococcus_lucimarinus","panicum_hallii","panicum_hallii_fil2","papaver_somniferum","phaseolus_vulgaris","physcomitrium_patens","pistacia_vera","populus_trichocarpa","prunus_avium","prunus_dulcis","prunus_persica","quercus_lobata","rosa_chinensis","saccharum_spontaneum","selaginella_moellendorffii","sesamum_indicum","setaria_italica","setaria_viridis","solanum_lycopersicum","solanum_tuberosum","solanum_tuberosum_rh8903916","sorghum_bicolor","theobroma_cacao","theobroma_cacao_criollo","trifolium_pratense","triticum_aestivum","triticum_aestivum_arinalrfor","triticum_aestivum_cadenza","triticum_aestivum_claire","triticum_aestivum_jagger","triticum_aestivum_julius","triticum_aestivum_lancer","triticum_aestivum_landmark","triticum_aestivum_mace","triticum_aestivum_mattis","triticum_aestivum_norin61","triticum_aestivum_paragon","triticum_aestivum_robigus","triticum_aestivum_stanley","triticum_aestivum_weebil","triticum_dicoccoides","triticum_spelta","triticum_turgidum","triticum_urartu","vigna_angularis","vigna_radiata","vitis_vinifera","zea_mays","ashbya_gossypii","aspergillus_clavatus","aspergillus_flavus","aspergillus_fumigatus","aspergillus_fumigatusa1163","aspergillus_nidulans","aspergillus_niger","aspergillus_oryzae","aspergillus_terreus","beauveria_bassiana","blumeria_graminis","botrytis_cinerea","colletotrichum_gloeosporioides","colletotrichum_graminicola","colletotrichum_higginsianum","colletotrichum_orbiculare","cryptococcus_neoformans","dothistroma_septosporum","fungi_ascomycota1_collection","fungi_ascomycota2_collection","fungi_ascomycota3_collection","fungi_ascomycota4_collection","fungi_basidiomycota1_collection","fungi_blastocladiomycota1_collection","fungi_chytridiomycota1_collection","fungi_entomophthoromycota1_collection","fungi_microsporidia1_collection","fungi_mucoromycota1_collection","fungi_neocallimastigomycota1_collection","fungi_rozellomycota1_collection","fusarium_culmorum","fusarium_fujikuroi","fusarium_graminearum","fusarium_oxysporum","fusarium_pseudograminearum","fusarium_solani","fusarium_verticillioides","gaeumannomyces_graminis","komagataella_pastoris","leptosphaeria_maculans","magnaporthe_oryzae","magnaporthe_poae","melampsora_laricipopulina","microbotryum_violaceum","neosartorya_fischeri","neurospora_crassa","phaeosphaeria_nodorum","puccinia_graminis","puccinia_graminisug99","puccinia_striiformis","puccinia_triticina","pyrenophora_teres","pyrenophora_triticirepentis","schizosaccharomyces_cryophilus","schizosaccharomyces_japonicus","schizosaccharomyces_octosporus","schizosaccharomyces_pombe","sclerotinia_sclerotiorum","sporisorium_reilianum","trichoderma_reesei","trichoderma_virens","tuber_melanosporum","ustilago_maydis","verticillium_dahliae","verticillium_dahliaejr2","yarrowia_lipolytica","zymoseptoria_tritici","actinia_equina_gca011057435","acyrthosiphon_pisum","adineta_vaga","aedes_aegypti_lvpagwg","aedes_albopictus","amphimedon_queenslandica","anopheles_albimanus","anopheles_arabiensis","anopheles_atroparvus","anopheles_christyi","anopheles_coluzzii","anopheles_coluzzii_ngousso","anopheles_culicifacies","anopheles_darlingi","anopheles_dirus","anopheles_epiroticus","anopheles_farauti","anopheles_funestus","anopheles_gambiae","anopheles_maculatus","anopheles_melas","anopheles_merus","anopheles_minimus","anopheles_quadriannulatus","anopheles_sinensis","anopheles_sinensis_china","anopheles_stephensi","anopheles_stephensi_indian","anoplophora_glabripennis","apis_mellifera","atta_cephalotes","belgica_antarctica","bemisia_tabaci_asiaii5","bemisia_tabaci_ssa1nig","bemisia_tabaci_ssa1ug","bemisia_tabaci_ssa2nig","bemisia_tabaci_ssa3nig","bemisia_tabaci_sweetpotug","biomphalaria_glabrata","bombus_impatiens","bombus_terrestris","bombyx_mori","branchiostoma_lanceolatum","brugia_malayi","caenorhabditis_brenneri","caenorhabditis_briggsae","caenorhabditis_japonica","caenorhabditis_remanei","capitella_teleta","cimex_lectularius","clytia_hemisphaerica_gca902728285","crassostrea_gigas","culex_quinquefasciatus","culicoides_sonorensis","danaus_plexippus","daphnia_magna","daphnia_pulex","dendroctonus_ponderosae","dinothrombium_tinctorium","drosophila_ananassae","drosophila_erecta","drosophila_grimshawi","drosophila_mojavensis","drosophila_persimilis","drosophila_pseudoobscura","drosophila_sechellia","drosophila_simulans","drosophila_virilis","drosophila_willistoni","drosophila_yakuba","folsomia_candida","glossina_austeni","glossina_brevipalpis","glossina_fuscipes","glossina_morsitans","glossina_pallidipes","glossina_palpalis","heliconius_melpomene","helobdella_robusta","hofstenia_miamia","ixodes_scapularis","ixodes_scapularis_ise6","lepeophtheirus_salmonis","leptotrombidium_deliense","lingula_anatina","loa_loa","lottia_gigantea","lucilia_cuprina","lutzomyia_longipalpis","mayetiola_destructor","megaselia_scalaris","melitaea_cinxia","mnemiopsis_leidyi","musca_domestica","nasonia_vitripennis","nematostella_vectensis","octopus_bimaculoides","onchocerca_volvulus","orchesella_cincta","pediculus_humanus","phlebotomus_papatasi","pristionchus_pacificus","rhodnius_prolixus","sarcoptes_scabiei","schistosoma_mansoni","solenopsis_invicta","stegodyphus_mimosarum","stomoxys_calcitrans","strigamia_maritima","strongylocentrotus_purpuratus","strongyloides_ratti","teleopsis_dalmanni","tetranychus_urticae","thelohanellus_kitauei","tigriopus_californicus_gca007210705","trialeurodes_vaporariorum_gca011764245","tribolium_castaneum","trichinella_spiralis","trichoplax_adhaerens","varroa_destructor_gca002443255","zootermopsis_nevadensis",	]
return list_genome_ensembl
}