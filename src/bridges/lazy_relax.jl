"""
Wraps an optimizer in a hierarchy which enables one `bridge`.
"""
function lazy_relax(optimizer, bridge::Type{<:MOIB.Constraint.AbstractBridge})
        function construct()
                model = Model(optimizer)
                JuMP.add_bridge(model, bridge)
                backend(model).optimizer
        end
        construct
end
