function fig = K_simulation(k,present_response,null_response)
fig =figure;
    for w = 1:length(k)
        for d = 1:length(present_response)
            
            % sort and calculate rank order of each patch response
            [rank_ord_response,I] = sort([present_response(d) null_response(d,:)],'descend');
            B = sort(unique([present_response(d) null_response(d)]),'descend');

            for i = 1:length(B)
                if i == 1
                    rank_order = ones(sum(rank_ord_response == B(i)),1);
                else
                    rank_current = ones(sum(rank_ord_response == B(i)),1).* (length(rank_order)+1);
                    rank_order = [rank_order; rank_current];
                end

            end
            
            % calculate probability
            present_rank = I(1);
            rank_order_p_raw = (length(rank_ord_response)-rank_order)/length(rank_ord_response);
            rank_order_p = (rank_order_p_raw.^k(w))/sum(rank_order_p_raw.^k(w));
            
            % plot probability

            subplot(length(present_response),length(k),(d-1).*length(k) + w),plot(0:4, rank_order_p,'.-b', 'MarkerSize',8,'LineWidth',1.5);
            hold on
            plot(present_rank-1, rank_order_p(present_rank),'d','MarkerFaceColor','red','MarkerEdgeColor','red');
            hold off

            set(gca,'XTick',[0:1:4],'XTickLabel',{num2str(rank_ord_response(1)),num2str(rank_ord_response(2)),num2str(rank_ord_response(3)),...
                num2str(rank_ord_response(4)),num2str(rank_ord_response(5))},'YTick',[0:0.2:0.6],'FontName','Arial','FontSize',12);
            ylim([0 0.6]);
            
            if w == 1
                ylabel('Probability');
            elseif w == 3
                xlabel('Decision x confidence');
            end
            
             if d == 1
                title(['k = ',num2str(k(w))]);
             end

        end
    end
end
