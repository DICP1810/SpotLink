VALUE_MAX_SCAN = 2000000
VALUE_APP_SCAN = 100000
VALUE_ILLEGAL = -716
VALUE_MAX_VALUE = 2000000

MS2_OUTPUT_TYPE = {
    'RAW': 0,
    'MS2_PROCESSED': 1
}

PSM_TYPE = {
    'RAW_PSM': 'rpsm',
    'FINE_SCORING_PSM':'fpsm',
    'SVM_PSM':'cxsvmpsm',
    'SITE_SCORING_PSM':'spsm',
    'SITE_FDR_PSM':'sfdrpsm',
}

PSM_OUTPUT_TYPE = {
    'ALL': 0,
    'TOP1_UNFILTERED': 1,
    'TOP1_FILTERED': 2,
}

GLOBAL_DIGEST_NAME = {
    'SPECIFIC': 0,
    'SEMI_SPECIFIC': 1,
    'UNSPECIFIC': 2
}

TOLERANCE_TYPE = {
    'ABSOLUTE_TOLERANCE': 0,
    'RELATIVE_TOLERANCE': 1
}

EXPIRATION_TIME = {'Year': 2419, 'Month': 12, 'Day': 31}

INFO_TO_USER_Staff = (
    '\n[SpotLink] Copyright \u00A9 2021. All rights reserved.',
    '\n[SpotLink] SpotLink is expired! Please send e-mail for the new version.',
    '\n[SpotLink] Warning! The current license will expired in 7 days. Please send e-mail for the new version.',
    '\n[SpotLink] Writing config file in the folder...',
    '\n[SpotLink] Finished!',
    '\n[SpotLink] Writing config file...',
    '\n[SpotLink] Run SpotLink with this command: python SpotLink [parameter file]',
    '\n[SpotLink] SpotLink whole workflow.',
    '\n[SpotLink] SpotLink fine scoring and PSM-FDR workflow.',
    '\n[SpotLink] SpotLink site scoring and site-FDR workflow.',)

CFLOW3_INFORMATION = (
    '[SpotLink] Start to read ini and ms2 file.',
    '[SpotLink] Start to generate protein & peptide index.',
    '[SpotLink] Digestion.',
    '[SpotLink] Generating Theoretical Spectra.',
    '[SpotLink] Saving indexes.',
    '[SpotLink] Loading exist indexes.',
    '[SpotLink] Start searching.',
    '[SpotLink] End searching.',
    '[SpotLink] Start Quality Control.',
    '[SpotLink] Start SVM classification on PSM.',
    '[SpotLink] Start Site Scoring.',
    '[SpotLink] Start sFDR calculation.'
)

CFLOW4_INFORMATION = (
    '[SpotLink] Start to read ini and ms2 file.',
    '[SpotLink] Start Quality Control.',
    '[SpotLink] Start SVM classification on PSM.',
    '[SpotLink] Start Site Scoring.',
    '[SpotLink] Start sFDR calculation.'
)

CFLOW5_INFORMATION = (
    '[SpotLink] Start to read ini and ms2 file.',
    '[SpotLink] Start Quality Control.',
    '[SpotLink] Start Site Scoring.',
    '[SpotLink] Start sFDR calculation.'
)

PRIME_NUMBER = (
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
    233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
    353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
    419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
    467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
    547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
    607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
    739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
    811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
    877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
    947, 953, 967, 971, 977, 983, 991, 997)

LOG_PRIME_NUMBER = (0.6931471805599453,
                    1.0986122886681098,
                    1.6094379124341003,
                    1.9459101490553132,
                    2.3978952727983707,
                    2.5649493574615367,
                    2.833213344056216,
                    2.9444389791664403,
                    3.1354942159291497,
                    3.367295829986474,
                    3.4339872044851463,
                    3.6109179126442243,
                    3.713572066704308,
                    3.7612001156935624,
                    3.8501476017100584,
                    3.970291913552122,
                    4.07753744390572,
                    4.110873864173311,
                    4.204692619390966,
                    4.2626798770413155,
                    4.290459441148391,
                    4.3694478524670215,
                    4.418840607796598,
                    4.48863636973214,
                    4.574710978503383,
                    4.61512051684126,
                    4.634728988229636,
                    4.672828834461906,
                    4.6913478822291435,
                    4.727387818712341,
                    4.844187086458591,
                    4.875197323201151,
                    4.919980925828125,
                    4.9344739331306915,
                    5.003946305945459,
                    5.017279836814924,
                    5.056245805348308,
                    5.093750200806762,
                    5.117993812416755,
                    5.153291594497779,
                    5.187385805840755,
                    5.198497031265826,
                    5.25227342804663,
                    5.262690188904886,
                    5.2832037287379885,
                    5.293304824724492,
                    5.351858133476067,
                    5.407171771460119,
                    5.424950017481403,
                    5.43372200355424,
                    5.4510384535657,
                    5.476463551931511,
                    5.484796933490655,
                    5.5254529391317835,
                    5.54907608489522,
                    5.572154032177765,
                    5.594711379601839,
                    5.602118820879701,
                    5.6240175061873385,
                    5.638354669333745,
                    5.645446897643238,
                    5.680172609017068,
                    5.726847747587197,
                    5.739792912179234,
                    5.746203190540153,
                    5.75890177387728,
                    5.802118375377063,
                    5.820082930352362,
                    5.849324779946859,
                    5.855071922202427,
                    5.8664680569332965,
                    5.883322388488279,
                    5.905361848054571,
                    5.921578419643816,
                    5.937536205082426,
                    5.948034989180646,
                    5.963579343618446,
                    5.983936280687191,
                    5.993961427306569,
                    6.013715156042802,
                    6.037870919922137,
                    6.042632833682381,
                    6.066108090103747,
                    6.07073772800249,
                    6.0844994130751715,
                    6.093569770045136,
                    6.1070228877422545,
                    6.124683390894205,
                    6.133398042996649,
                    6.137727054086234,
                    6.1463292576688975,
                    6.171700597410915,
                    6.18826412308259,
                    6.19644412779452,
                    6.212606095751519,
                    6.220590170099739,
                    6.2324480165505225,
                    6.255750041753367,
                    6.259581464064923,
                    6.293419278846481,
                    6.304448802421981,
                    6.322565239927284,
                    6.333279628139691,
                    6.343880434126331,
                    6.3473892096560105,
                    6.3578422665081,
                    6.375024819828097,
                    6.385194398997726,
                    6.395261598115449,
                    6.398594934535208,
                    6.408528791059498,
                    6.418364935936212,
                    6.424869023905388,
                    6.428105272684596,
                    6.447305862541213,
                    6.46302945692067,
                    6.466144724237619,
                    6.472346294500901,
                    6.481577129276431,
                    6.490723534502507,
                    6.493753839851686,
                    6.511745329644728,
                    6.517671272912275,
                    6.52649485957079,
                    6.53813982376767,
                    6.55250788703459,
                    6.5638555265321274,
                    6.577861357721047,
                    6.588926477533519,
                    6.597145701886651,
                    6.6052979209482015,
                    6.610696044717759,
                    6.621405651764134,
                    6.6293632534374485,
                    6.634633357861686,
                    6.645090969505644,
                    6.650279048587422,
                    6.668228248417403,
                    6.680854678790215,
                    6.695798917058491,
                    6.698268054115413,
                    6.710523109452428,
                    6.71295620067707,
                    6.717804695023691,
                    6.720220155135295,
                    6.732210706467206,
                    6.748759547491679,
                    6.75343791859778,
                    6.755768921984255,
                    6.760414691083428,
                    6.776506992372183,
                    6.78105762593618,
                    6.78332520060396,
                    6.787844982309579,
                    6.810142450115136,
                    6.814542897259958,
                    6.823286122355687,
                    6.834108738813838,
                    6.842683282238422,
                    6.846943139585379,
                    6.853299093186078,
                    6.859614903654202,
                    6.874198495453294,
                    6.878326468291325,
                    6.884486652042782,
                    6.890609120147166,
                    6.898714534329988,
                    6.904750769961838)
