#####

# Import pyidlrpc.
import pyidlrpc as idl

# Execute a command in IDL.
idl.execute( $
maplist = HP_AKARI_list.txt
inputlist = AME_galcoord.txt
radius = 60
rinner = 1.20
router = 2.0
;;;See comments about how to properly prepare non-Planck data to be used by this routine...
multiepoch_photometry, inputlist, maplist=maplist, radius=radius, rinner=rinner, router=router, /galactic, /noise_model
END $
)
