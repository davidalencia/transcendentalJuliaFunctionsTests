using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
  "fun"
      required = true
  "--rng_start", "-s"
      required = true
      arg_type = Float64
  "--rng_end", "-e"
      required = true
      arg_type = Float64
  "--epochs"
      help = "another option with an argument"
      arg_type = Int
      default = 1000
end

args = parse_args(ARGS, s)

println(ARGS)