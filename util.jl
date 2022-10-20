function get_debug()::Bool
    for arg in ARGS
        if arg == "debug" return true end
    end

    return false
end
