import { MolstarLoadingExtension } from 'molstar/lib/extensions/mvs/load';

export const IsHiddenCustomStateExtension: MolstarLoadingExtension<{}> = {
    id: 'is-hidden-custom-state',
    description: 'Allow updating initial visibility of nodes',
    createExtensionContext: () => ({}),
    action: (updateTarget, node) => {
        if (!node.custom || !node.custom?.is_hidden) return;
        updateTarget.update.to(updateTarget.selector).updateState({ isHidden: true });
    },
};