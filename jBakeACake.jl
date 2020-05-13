module jBakeACake

export NeedIngredients, Gather, Add, MixBatter, Bake
export Ingredient, eggs, flour, sugar, water

abstract type State end
abstract type Transition end

invalidState = ArgumentError
@enum Ingredient eggs flour sugar water butter baking_soda

struct NeedIngredients <: State
  gathered::Set{Ingredient}
  NeedIngredients() = new(Set())
end

struct Gather <: Transition
  ingredient::Ingredient
end

struct MixBatter <: State
  remaining::Array{Ingredient}
  MixBatter() = new( [flour, sugar, baking_soda, water, eggs, butter] )
end

struct Add <: Transition
  ingredient::Ingredient
end

struct RuinedBatter <: State end
struct RawBatter <: State end
struct TastyCake <: State end
struct BurntCake <: State end

struct Bake <: Transition
  temp::UInt64
end

step(::State, ::Transition) = throw(invalidState)

function step(state::NeedIngredients, transition::Gather)
  push!(state.gathered, transition.ingredient)
  state.gathered != Set{instances(Ingredient)} ? state : MixBatter()
end

function step(state::MixBatter, transition::Add)
  if state.remaining[1] == transition.ingredient
    popfirst!(state.remaining)
    isempty(state.remaining) ? RawBatter() : state
  else
    println("Ingredients added in wrong order")
    RuinedBatter()
  end
end

function Bake(::RawBatter, transition::Bake) transition.temp > 350 ? BurntCake() : TastyCake() end

end  # module jBakeACake
