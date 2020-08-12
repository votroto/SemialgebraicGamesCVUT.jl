function available_result(backend)
	pstat = MOI.get(backend, MOI.PrimalStatus())
	dstat = MOI.get(backend, MOI.DualStatus())
	return pstat != MOI.UNKNOWN_RESULT_STATUS &&
	       pstat != MOI.NO_SOLUTION &&
	       dstat != MOI.UNKNOWN_RESULT_STATUS &&
	       dstat != MOI.NO_SOLUTION
end