;;;;;; stepH_haper.pro

maplist = HP_AKARI_list.txt
inputlist = AME_galcoord.txt
radius = 60
galactic
rinner = 1.20
router = 2.0

multiepoch_photometry, inputlist, maplist=maplist, radius=radius, galactic=galactic, decimal=decimal, rinner=rinner, router=router
END
