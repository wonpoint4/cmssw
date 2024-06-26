#ifndef FWCore_Framework_one_EDFilter_h
#define FWCore_Framework_one_EDFilter_h
// -*- C++ -*-
//
// Package:     FWCore/Framework
// Class  :     edm::one::EDFilter
//
/**\class edm::one::EDFilter EDFilter.h "FWCore/Framework/interface/one/EDFilter.h"

 Description: [one line class summary]

 Usage:
    <usage>

*/
//
// Original Author:  Chris Jones
//         Created:  Thu, 09 May 2013 19:53:55 GMT
//

// system include files

// user include files
#include "FWCore/Framework/interface/one/filterAbilityToImplementor.h"

// forward declarations
namespace edm {
  namespace one {
    template <typename... T>
    class EDFilter : public filter::AbilityToImplementor<T>::Type..., public virtual EDFilterBase {
    public:
      static_assert(not(CheckAbility<module::Abilities::kRunCache, T...>::kHasIt and
                        CheckAbility<module::Abilities::kOneWatchRuns, T...>::kHasIt),
                    "Cannot use both WatchRuns and RunCache");
      static_assert(not(CheckAbility<module::Abilities::kLuminosityBlockCache, T...>::kHasIt and
                        CheckAbility<module::Abilities::kOneWatchLuminosityBlocks, T...>::kHasIt),
                    "Cannot use both WatchLuminosityBlocks and LuminosityBLockCache");
      EDFilter() = default;
      EDFilter(const EDFilter&) = delete;
      const EDFilter& operator=(const EDFilter&) = delete;

      //virtual ~EDFilter();

      // ---------- const member functions ---------------------
      bool wantsProcessBlocks() const noexcept final { return WantsProcessBlockTransitions<T...>::value; }
      bool wantsInputProcessBlocks() const noexcept final { return WantsInputProcessBlockTransitions<T...>::value; }
      bool wantsGlobalRuns() const noexcept final { return WantsGlobalRunTransitions<T...>::value; }
      bool wantsGlobalLuminosityBlocks() const noexcept final {
        return WantsGlobalLuminosityBlockTransitions<T...>::value;
      }

      bool hasAbilityToProduceInBeginProcessBlocks() const final {
        return HasAbilityToProduceInBeginProcessBlocks<T...>::value;
      }
      bool hasAbilityToProduceInEndProcessBlocks() const final {
        return HasAbilityToProduceInEndProcessBlocks<T...>::value;
      }

      bool hasAbilityToProduceInBeginRuns() const final { return HasAbilityToProduceInBeginRuns<T...>::value; }
      bool hasAbilityToProduceInEndRuns() const final { return HasAbilityToProduceInEndRuns<T...>::value; }

      bool hasAbilityToProduceInBeginLumis() const final { return HasAbilityToProduceInBeginLumis<T...>::value; }
      bool hasAbilityToProduceInEndLumis() const final { return HasAbilityToProduceInEndLumis<T...>::value; }

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      SerialTaskQueue* globalRunsQueue() final { return globalRunsQueue_.queue(); }
      SerialTaskQueue* globalLuminosityBlocksQueue() final { return globalLuminosityBlocksQueue_.queue(); }

    private:
      // ---------- member data --------------------------------
      impl::OptionalSerialTaskQueueHolder<WantsSerialGlobalRunTransitions<T...>::value> globalRunsQueue_;
      impl::OptionalSerialTaskQueueHolder<WantsSerialGlobalLuminosityBlockTransitions<T...>::value>
          globalLuminosityBlocksQueue_;
    };

  }  // namespace one
}  // namespace edm

#endif
