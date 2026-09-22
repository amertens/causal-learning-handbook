# Citation additions — verification checklist

34 foundational references were added to `book.bib` on 2026-07-12 by the citation pass.

**VERIFIED 2026-07-12 against Crossref (api.crossref.org).** 30 of 34 entries were confirmed by DOI: volume, issue, and page ranges match what is in `book.bib` (a few entries returned only a start page or no pages from Crossref but had volume/issue/year confirmed; their full page ranges are from standard citations). The remaining 4 have no reliable Crossref DOI and were confirmed against standard citations: Manski 1990 (Am Econ Rev 80(2):319-323), Aalen & Johansen 1978 (Scand J Statist 5(3):141-150), Pearl 2001 (Proc. UAI 17, pp. 411-420), and the `benkeser2017survtmle` software `@manual` (no pages). **No corrections were needed.** The "Verify" column below is retained for reference; every row has been checked.

| Key | Reference | Verify |
|-----|-----------|--------|
| robins1986gformula | Robins (1986), Mathematical Modelling 7(9-12):1393-1512 | vol/pages |
| robins2000msm | Robins, Hernán & Brumback (2000), Epidemiology 11(5):550-560 | ok, spot-check pages |
| cole2008ipw | Cole & Hernán (2008), Am J Epidemiol 168(6):656-664 | pages |
| greenland1999collapsibility | Greenland, Robins & Pearl (1999), Statist Sci 14(1):29-46 | pages |
| suissa2008immortal | Suissa (2008), Am J Epidemiol 167(4):492-499 | pages |
| vanderweele2017evalue | VanderWeele & Ding (2017), Ann Intern Med 167(4):268-274 | pages |
| ding2016sensitivity | Ding & VanderWeele (2016), Epidemiology 27(3):368-377 | pages |
| rothman1976synergy | Rothman (1976), Am J Epidemiol 103(5):506-511 | pages |
| athey2019grf | Athey, Tibshirani & Wager (2019), Ann Statist 47(2):1148-1178 | pages |
| kunzel2019metalearners | Künzel, Sekhon, Bickel & Yu (2019), PNAS 116(10):4156-4165 | pages |
| robinson1988 | Robinson (1988), Econometrica 56(4):931-954 | pages |
| wolpert1992stacking | Wolpert (1992), Neural Networks 5(2):241-259 | pages |
| breiman1996stacking | Breiman (1996), Machine Learning 24(1):49-64 | pages |
| dudoit2005asymptotics | Dudoit & van der Laan (2005), Statistical Methodology 2(2):131-154 | vol/pages |
| westreichcole2010positivity | Westreich & Cole (2010), Am J Epidemiol 171(6):674-677 | pages |
| crump2009trimming | Crump, Hotz, Imbens & Mitnik (2009), Biometrika 96(1):187-199 | pages |
| lipsitch2010negcontrol | Lipsitch, Tchetgen Tchetgen & Cohen (2010), Epidemiology 21(3):383-388 | pages |
| miao2018proximal | Miao, Geng & Tchetgen Tchetgen (2018), Biometrika 105(4):987-993 | pages |
| manski1990bounds | Manski (1990), Am Econ Rev 80(2):319-323 | pages |
| aalen1978 | Aalen & Johansen (1978), Scand J Statist 5(3):141-150 | pages |
| finegray1999 | Fine & Gray (1999), JASA 94(446):496-509 | pages |
| putter2007competing | Putter, Fiocco & Geskus (2007), Stat Med 26(11):2389-2430 | pages |
| benkeser2017survtmle | Benkeser & Hejazi, survtmle R package (2017) | version/year |
| stensrud2022separable | Stensrud, Young, Didelez, Robins & Hernán (2022), JASA 117(537):175-183 | pages |
| robinsgreenland1992 | Robins & Greenland (1992), Epidemiology 3(2):143-155 | pages |
| pearl2001mediation | Pearl (2001), Proc. UAI 17:411-420 | pages |
| vanderweele2014mediation | VanderWeele, Vansteelandt & Robins (2014), Epidemiology 25(2):300-306 | pages |
| tchetgen2012mediation | Tchetgen Tchetgen & Shpitser (2012), Ann Statist 40(3):1816-1845 | pages |
| imbens1994late | Imbens & Angrist (1994), Econometrica 62(2):467-475 | pages |
| colnet2024combining | Colnet et al. (2024), Statist Sci 39(1):165-191 | vol/pages |
| dahabreh2020transport | Dahabreh et al. (2020), Stat Med 39(14):1999-2014 | pages |
| bareinboim2016fusion | Bareinboim & Pearl (2016), PNAS 113(27):7345-7352 | pages |
| rudolph2017transport | Rudolph & van der Laan (2017), JRSS-B 79(5):1509-1525 | pages |
| hahn1998efficiency | Hahn (1998), Econometrica 66(2):315-331 | pages |

Quick verification: `doi.org` or the journal's site for each; or re-run `bibtex_utils.py fetch --doi <doi>` when the network allows.
