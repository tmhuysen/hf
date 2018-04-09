#ifndef HF_SCFSOLVERTYPE_HPP
#define HF_SCFSOLVERTYPE_HPP

namespace hf{
namespace solver{


/**
 *  An enum class for the implemented SCF solver types.
 */
enum class SCFSolverType {
    PLAIN,
    DIIS
};


}  // namespace solver
}  // namespace hf



#endif //HF_SCFSOLVERTYPE_HPP
