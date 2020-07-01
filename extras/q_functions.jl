function huber_loss(x, δ::Float64=1.0)
    if abs(x) < δ
        return 0.5*x^2
    else
        return δ*(abs(x) - 0.5*δ)
    end
end

function POMDPs.discount(pomdp::PowerGridEnv)
    return 0.8
end

function q_learning_loss(solver, env, active_q, target_q, s_batch, a_batch, r_batch, sp_batch, done_batch, importance_weights)
    q_values = active_q(s_batch)
    q_sa = diag(view(q_values, a_batch, :))
    if solver.double_q
        target_q_values = target_q(sp_batch)
        qp_values = active_q(sp_batch)
        best_a = Flux.onecold(qp_values)
        q_sp_max = diag(view(target_q_values, best_a, :))
    else
        q_sp_max = @view maximum(target_q(sp_batch), dims=1)[:]
    end
    q_targets = r_batch .+ convert(Vector{Float32}, (1.0 .- done_batch) .* discount(env.problem)).*q_sp_max
    td_tracked = q_sa .- q_targets
    loss_tracked = mean(huber_loss, importance_weights.*td_tracked)
    return loss_tracked, td_tracked
end
