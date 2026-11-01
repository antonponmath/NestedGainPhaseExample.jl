module Timing

export timed

_MSG = Dict(
    :PREPARE => "Preliminary calculations",
    :SIMPLE => " μ₁G₁ ↺ (I + Δ) with ‖Δ‖ ≤ μ₂‖G₂‖",
    :PARALLEL => " μ₁G₁ ↺ (I + μ₂G₂)",
    :INNER_UNITY => "(μ₁G₁ ↺ I) ↺ μ₂G₂",
    :OUTER_UNITY => "(μ₁G₁ ↺ μ₂G₂) ↺ I",
    :SAMPLING => "Parameter sampling",
    :DRAWING => "Drawing",
)

_TMP = "...working..."

_msg_finish = T -> "\b"^length(_TMP) * "done in $(round(T, digits = 3)) seconds."

_max_len = maximum([length(msg) for (_, msg) in _MSG])

for (key, msg) in _MSG
    _MSG[key] = rpad(msg, _max_len + 4) * _TMP
end

function timed(fn, key::Symbol)
    print(_MSG[key])
    t0 = time()
    result = fn()
    println(_msg_finish(time() - t0))
    return result
end

end # module Timing
